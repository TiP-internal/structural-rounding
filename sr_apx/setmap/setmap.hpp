
#ifndef SETMAP_H
#define SETMAP_H

// needed for std::ceil
#include <cmath>

// needed for std::max
#include <algorithm>

// needed for std::out_of_range
#include <stdexcept>

namespace sr_apx {

/**
	The setmap_detail namespace contains implementation specifics of set and map.
	The classes/functions defined within should not be used directly.
*/
namespace setmap_detail {

/**
	The HashWrap class wraps a given hash function so that
	it can be called with both keys and (key,value) pairs.

	ValueType is presumed to be std::pair<KeyType, T>.

	Hash must be a function object which defines operator(),
	accepting a KeyType and producing a std::size_t in the manner
	of std::hash<KeyType>.
*/
template<class KeyType, class ValueType, class Hash>
struct HashWrap: private Hash {
	// HashWrap inherits from Hash to use empty base optimization

	HashWrap() = default;
	HashWrap(const Hash& hash): Hash(hash) {}

	std::size_t operator()(const KeyType& key) const {
		return Hash::operator()(key);
	}

	std::size_t operator()(const ValueType& value) const {
		return Hash::operator()(value.first);
	}
};

/**
	The EqualWrap class wraps a given equals function so that
	it can be called with keys or (key,value) pairs for both
	arguments.

	ValueType is presumed to be std::pair<KeyType, T>.

	Equal must be a function object which defines operator(),
	accepting two KeyType's and producing a bool in the manner
	of std::equal_to<KeyType>.
*/
template<class KeyType, class ValueType, class Equal>
struct EqualWrap: private Equal {
	// EqualWrap inherits from Equal to use empty base optimization

	EqualWrap() = default;
	EqualWrap(const Equal& equal): Equal(equal) {}

	bool operator()(const KeyType& lhs, const KeyType& rhs) const {
		return Equal::operator()(lhs, rhs);
	}

	bool operator()(const KeyType& lhs, const ValueType& rhs) const {
		return Equal::operator()(lhs, rhs.first);
	}

	bool operator()(const ValueType& lhs, const KeyType& rhs) const {
		return Equal::operator()(lhs.first, rhs);
	}

	bool operator()(const ValueType& lhs, const ValueType& rhs) const {
		return Equal::operator()(lhs.first, rhs.first);
	}
};

/**
	The Bucket class combines storage for a ValueType object
	with metadata needed for robin hood hashing.

	The value field is placed in an anonymous union to prevent
	it from being default initialized when a bucket is created.
	The value field must be manually constructed using set_value.
	Destruction is handled automatically.
*/
template<class ValueType>
struct Bucket {
	// metadata field for robin hood hashing algorithm
	// negative age signifies empty
	int age = -1;

	// storage space for items
	// placed in a union to prevent default construction
	union {
		ValueType value;
	};

	Bucket() noexcept {}

	~Bucket() noexcept {
		clear();
	}

	void clear() noexcept {
		// only calls the destructor if not empty
		if (age != -1) {
			value.~ValueType();
			age = -1;
		}
	}

	bool empty() const noexcept {
		return age == -1;
	}

	template<class... Args>
	void set_value(Args&&... args) {
		// placement new, constructs a new object at the specified address
		new(&value) ValueType(std::forward<Args>(args)...);
	}
};

/**
	The table class is the shared base class of both set and map.
	It implements a hash table using open addressing with robin hood
	hashing for collision resolution.

	ValueType is presumed to be the same as KeyType or
	std::pair<KeyType, T>.

	Hash must be a function object which defines operator().
	It should return a std::size_t and accept either a KeyType
	or ValueType.

	Equal must be a function object which defines operator().
	It should return a bool and accept 2 arguments, either KeyType's
	or ValueType's.
*/
template<
	class KeyType,
	class ValueType,
	class Hash,
	class Equal>
class table: private Hash, private Equal {
	// table inherits from Hash and Equal to use empty base optimization

public:
	// iterator type for sets and maps
	// declaring in advance for using declarations
	template<class IterType>
	class table_iterator;

	using key_type = KeyType;
	using value_type = ValueType;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;
	using hasher = Hash;
	using key_equal = Equal;
	using reference = value_type&;
	using const_reference = const value_type&;
	using pointer = value_type*;
	using const_pointer = const value_type*;
	using iterator = table_iterator<value_type>;
	using const_iterator = table_iterator<const value_type>;

private:
	using bucket = Bucket<value_type>;

	// maximum load factor of a set/map before it rehashes into a larger container
	static constexpr float max_lf = 0.5;

	// storage array
	bucket* buckets = nullptr;

	// number of elements
	size_type load = 0;

	// size of buckets array
	size_type capacity = 0;

public:

	/**
		The table_iterator class implements a foward iterator for tables.
		It iterates backwards over the allocated buckets.

		IterType is either value_type or const value_type, allowing one class
		to implement both iterators and const_iterators.
	*/
	template<class IterType>
	class table_iterator {
		// reference to the storage array
		bucket* base = nullptr;

	public:
		// current index in the storage array
		size_type index = -1;

		using iterator_category = std::forward_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type*;
		using reference = value_type&;

		table_iterator(bucket* b, size_type i): base(b), index(i) {}

		friend bool operator==(const table_iterator& lhs, const table_iterator& rhs) {
			return lhs.index == rhs.index;
		}

		friend bool operator!=(const table_iterator& lhs, const table_iterator& rhs) {
			return !(lhs == rhs);
		}

		table_iterator& operator++() {
			if (index == -1) {
				return *this;
			}

			do {
				--index;
			}
			while (index != -1 && base[index].empty());

			return *this;
		}

		table_iterator operator++(int) {
			table_iterator copy(*this);
			++*this;
			return copy;
		}

		IterType& operator*() const {
			return base[index].value;
		}

		IterType* operator->() const {
			return &(base[index].value);
		}

		operator table_iterator<const value_type>() const {
			return table_iterator<const value_type>(base, index);
		}
	};

	table(): Hash(), Equal() {}

	explicit table(size_type n, const Hash& hash = Hash(), const Equal& equal = Equal()): Hash(hash), Equal(equal) {
		if(n > 0) {
			// allocates space for n elements, rather than n buckets
			reserve(n);
		}
	}

	template<class InputIt>
	table(InputIt first, InputIt last, size_type n = 0, const Hash& hash = Hash(), const Equal& equal = Equal()): table(n, hash, equal) {
		insert(first, last);
	}

	table(std::initializer_list<value_type> init, size_type n = 0, const Hash& hash = Hash(), const Equal& equal = Equal()): table(n, hash, equal) {
		if (init.size() > n) {
			reserve(init.size());
		}
		insert(init.begin(), init.end());
	}

	table(const table& other): Hash(other), Equal(other) {
		if (other.load > 0) {
			rehash_impl(other.capacity);
			insert_nocheck(other.begin(), other.end());
		}
	}

	table(table&& other) noexcept: table() {
		swap(other);
	}

	~table() {
		if (buckets != nullptr) {
			delete[] buckets;
			buckets = nullptr;
			capacity = 0;
			load = 0;
		}
	}

	table& operator=(const table& other) {
		if (this == &other) {
			return *this;
		}

		Hash::operator=(other);
		Equal::operator=(other);
		clear();

		if (other.load > 0) {
			rehash_impl(other.capacity);
			insert_nocheck(other.begin(), other.end());
		}

		return *this;
	}

	table& operator=(table&& other) noexcept {
		swap(other);
		return *this;
	}

	table& operator=(std::initializer_list<value_type> ilist) {
		clear();
		reserve(ilist.size());
		insert(ilist.begin(), ilist.end());
		return *this;
	}

	const Hash& hash_function() const {
		return static_cast<const Hash&>(*this);
	}

	const Equal& key_eq() const {
		return static_cast<const Equal&>(*this);
	}

	iterator begin() noexcept {
		return ++iterator(buckets, capacity);
	}

	const_iterator begin() const noexcept {
		return cbegin();
	}

	const_iterator cbegin() const noexcept {
		return ++iterator(buckets, capacity);
	}

	iterator end() noexcept {
		return iterator(nullptr, -1);
	}

	const_iterator end() const noexcept {
		return cend();
	}

	const_iterator cend() const noexcept {
		return const_iterator(nullptr, -1);
	}

	bool empty() const noexcept {
		return load == 0;
	}

	size_type size() const noexcept {
		return load;
	}

	size_type max_size() const noexcept {
		return std::numeric_limits<size_type>::max() / sizeof(bucket);
	}

	void clear() noexcept {
		if (buckets == nullptr) {
			return;
		}

		for (bucket* it = buckets; it < buckets + capacity; ++it) {
			it->clear();
		}

		load = 0;
	}

	std::pair<iterator, bool> insert(const value_type& value) {
		value_type temp(value);
		return insert_impl(std::move(temp));
	}

	std::pair<iterator, bool> insert(value_type&& value) {
		return insert_impl(std::move(value));
	}

	iterator insert(const_iterator hint, const value_type& value) {
		value_type temp(value);
		return insert_impl(std::move(temp)).first;
	}

	iterator insert(const_iterator hint, value_type&& value) {
		return insert_impl(std::move(value)).first;
	}

	template<class InputIt>
	void insert(InputIt first, InputIt last) {
		for (InputIt it = first; it != last; ++it) {
			value_type temp(*it);
			insert_impl(std::move(temp));
		}
	}

	void insert(std::initializer_list<value_type> ilist) {
		insert(ilist.begin(), ilist.end());
	}

	template<class... Args>
	std::pair<iterator, bool> emplace(Args&&... args) {
		value_type temp(std::forward<Args>(args)...);
		return insert_impl(std::move(temp));
	}

	template<class... Args>
	iterator emplace_hint(const_iterator hint, Args&&... args) {
		value_type temp(std::forward<Args>(args)...);
		return insert_impl(std::move(temp)).first;
	}

	iterator erase(const_iterator pos) {
		erase_impl(pos.index);
		return ++iterator(pos);
	}

	iterator erase(const_iterator first, const_iterator last) {
		// TODO: inefficient, shifts forward on every erase call
		for (const_iterator it = first; it != last; ++it) {
			erase_impl(it.index);
		}

		return ++iterator(last);
	}

	size_type erase(const key_type& key) {
		return erase_impl(find_impl(key));
	}

	void swap(table& other) {
		std::swap(buckets, other.buckets);
		std::swap(load, other.load);
		std::swap(capacity, other.capacity);

		std::swap(static_cast<Hash&>(*this), static_cast<Hash&>(other));
		std::swap(static_cast<Equal&>(*this), static_cast<Equal&>(other));
	}

	friend void swap(table& lhs, table& rhs) {
		lhs.swap(rhs);
	}

	size_type count(const key_type& key) const {
		size_type index = find_impl(key);
		return index == -1 ? 0 : 1;
	}

	iterator find(const key_type& key) {
		size_type index = find_impl(key);
		return index == -1 ? end() : iterator(buckets, index);
	}

	const_iterator find(const key_type& key) const {
		size_type index = find_impl(key);
		return index == -1 ? cend() : const_iterator(buckets, index);
	}

	bool contains(const key_type& key) const {
		return find_impl(key) != -1;
	}

	std::pair<iterator, iterator> equal_range(const key_type& key) {
		size_type index = find_impl(key);
		if (index != -1) {
			return {iterator(buckets, index), ++iterator(buckets, index)};
		}

		return {end(), end()};
	}

	std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const {
		size_type index = find_impl(key);
		if (index != -1) {
			return {const_iterator(buckets, index), ++const_iterator(buckets, index)};
		}

		return {cend(), cend()};
	}

	size_type bucket_count() const {
		return capacity;
	}

	size_type max_bucket_count() const {
		return std::numeric_limits<size_type>::max() / sizeof(bucket);
	}

	float load_factor() const {
		return load / capacity;
	}

	float max_load_factor() const {
		return max_lf;
	}

	void max_load_factor(float ml) {}

	/**
		Rehashes the container so that it contains at least count buckets.
		Ensures that the capacity is large enough for the load with the maximum load factor.
	*/
	void rehash(size_type count) {
		count = std::max(count, static_cast<size_type>(std::ceil(load / max_lf)));
		if (count <= capacity) {
			return;
		}

		--count;
		count |= count >> 1;
		count |= count >> 2;
		count |= count >> 4;
		count |= count >> 8;
		count |= count >> 16;
		count |= count >> 32;
		++count;

		rehash_impl(count);
	}

	/**
		Rehashes the container so that it can contain at least count elements.
	*/
	void reserve(size_type count) {
		rehash(std::ceil(count / max_lf));
	}

protected:
	void rehash_impl(size_type count) {
		count = std::max(count, static_cast<size_type>(64));

		if (count <= capacity) {
			return;
		}

		bucket* new_buckets = new bucket[count];

		std::swap(buckets, new_buckets);
		std::swap(capacity, count);
		load = 0;

		if (new_buckets == nullptr) {
			return;
		}

		for (bucket* it = new_buckets; it < new_buckets + count; ++it) {
			if (!it->empty()) {
				insert_impl_nocheck(Hash::operator()(it->value) & (capacity - 1), 0, std::move(it->value));
				++load;
				it->clear();
			}
		}

		delete[] new_buckets;
	}

	size_type erase_impl(size_type index) {
		if (index == -1 || buckets[index].empty()) {
			return 0;
		}

		size_type mask = capacity - 1;

		buckets[index].clear();
		--load;

		size_type prev_index = index;
		index = (index + 1) & mask;
		while (buckets[index].age > 0) {
			buckets[prev_index].set_value(std::move(buckets[index].value));
			buckets[prev_index].age = buckets[index].age - 1;

			buckets[index].age = -1;

			prev_index = index;
			index = (index + 1) & mask;
		}

		return 1;
	}

	size_type find_impl(const key_type& key) const {
		if (capacity == 0) {
			return -1;
		}

		size_type mask = capacity - 1;
		size_type index = Hash::operator()(key) & mask;
		int age = 0;

		while (age <= buckets[index].age) {
			if (Equal::operator()(key, buckets[index].value)) {
				return index;
			}

			index = (index + 1) & mask;
			++age;
		}

		return -1;
	}

	std::pair<iterator, bool> insert_impl(value_type&& value) {
		if (load + 1 > capacity * max_lf) {
			rehash_impl(capacity << 1);
		}

		size_type mask = capacity - 1;
		size_type index = Hash::operator()(value) & mask;

		int age = 0;

		while (age <= buckets[index].age) {
			if (Equal::operator()(value, buckets[index].value)) {
				return {iterator(buckets, index), false};
			}

			index = (index + 1) & mask;
			++age;
		}

		++load;

		if (buckets[index].empty()) {
			buckets[index].set_value(std::move(value));
			buckets[index].age = age;
		}
		else {
			std::swap(age, buckets[index].age);
			std::swap(value, buckets[index].value);
			insert_impl_nocheck((index + 1) & mask, age + 1, std::move(value));
		}

		return {iterator(buckets, index), true};
	}

	template<class InputIt>
	void insert_nocheck(InputIt first, InputIt last) {
		for (InputIt it = first; it != last; ++it) {
			value_type temp(*it);
			insert_impl_nocheck(Hash::operator()(temp) & (capacity - 1), 0, std::move(temp));
			++load;
		}
	}

	void insert_impl_nocheck(size_type index, int age, value_type&& value) {
		size_type mask = capacity - 1;
		while (!buckets[index].empty()) {
			if (age > buckets[index].age) {
				std::swap(age, buckets[index].age);
				std::swap(value, buckets[index].value);
			}

			index = (index + 1) & mask;
			++age;
		}

		buckets[index].set_value(std::move(value));
		buckets[index].age = age;
	}
};

}

/**
	The set class implements a hash set based on the table class above.

	Hash should be a function object which implements operator(), taking a Key
	and producing a std::size_t used as an index for storing the key.

	KeyEqual should be a function object with implements operator(), taking two
	Key objects and producing a bool which is true if the keys are equal.
*/
template<
	class Key,
	class Hash = std::hash<Key>,
	class KeyEqual = std::equal_to<Key>>
class set: public setmap_detail::table<Key, Key, Hash, KeyEqual> {
	using table = typename setmap_detail::table<Key, Key, Hash, KeyEqual>;

public:
	// gives set access to table constructors
	using table::table;

	friend bool operator==(const set& lhs, const set& rhs) {
		if (lhs.size() != rhs.size()) {
			return false;
		}

		for (typename table::const_iterator it = lhs.begin(); it != lhs.end(); ++it) {
			if (!rhs.contains(*it)) {
				return false;
			}
		}

		return true;
	}

	friend bool operator!=(const set& lhs, const set& rhs) {
		return !(lhs == rhs);
	}
};

/**
	The map class implements a hash map based on the table class above.
	The elements stored in the table are (key, value) pairs with the type
	std::pair<Key, T>.

	Hash should be a function object which implements operator(), taking a Key
	and producing a std::size_t used as an index for storing the key.
	The Hash object is wrapped so that it can also be called with (key, value)
	pairs.

	KeyEqual should be a function object with implements operator(), taking two
	Key objects and producing a bool which is true if the keys are equal.
	The KeyEqual object is wrapped so that it can also be called with (key, value)
	pairs.
*/
template<
	class Key,
	class T,
	class Hash = std::hash<Key>,
	class KeyEqual = std::equal_to<Key>>
class map: public setmap_detail::table<Key, std::pair<Key, T>, setmap_detail::HashWrap<Key, std::pair<Key, T>, Hash>, setmap_detail::EqualWrap<Key, std::pair<Key, T>, KeyEqual>> {
	using table = typename setmap_detail::table<Key, std::pair<Key, T>, setmap_detail::HashWrap<Key, std::pair<Key, T>, Hash>, setmap_detail::EqualWrap<Key, std::pair<Key, T>, KeyEqual>>;

public:
	// gives map access to table constructors
	using table::table;

	friend bool operator==(const map& lhs, const map& rhs) {
		if (lhs.size() != rhs.size()) {
			return false;
		}

		for (typename table::const_iterator lit = lhs.begin(); lit != lhs.end(); ++lit) {
			typename table::const_iterator rit = rhs.find(lit->first);
			if (rit == rhs.end() || lit->second != rit->second) {
				return false;
			}
		}

		return true;
	}

	friend bool operator!=(const map& lhs, const map& rhs) {
		return !(lhs == rhs);
	}

	T& operator[](const Key& key) {
		Key temp(key);
		return operator[](std::move(temp));
	}

	T& operator[](Key&& key) {
		std::pair<Key, T> val(std::move(key), std::move(T()));
		return table::insert_impl(std::move(val)).first->second;
	}

	T& at(const Key& key) {
		typename table::iterator it = table::find(key);
		if (it == table::end()) {
			throw std::out_of_range("key not in map");
		}

		return it->second;
	}

	const T& at(const Key& key) const {
		typename table::const_iterator it = table::find(key);
		if (it == table::end()) {
			throw std::out_of_range("key not in map");
		}

		return it->second;
	}

	// TODO insert_or_assign
	// TODO try_emplace
};

using Set = set<int>;
template<class T>
using Map = map<int, T>;

}

#endif
