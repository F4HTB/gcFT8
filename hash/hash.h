#ifndef _HASH_H_
#define _HASH_H_

#define hash_CAPACITY 50000

#ifdef __cplusplus
extern "C"
{
#endif

	typedef struct Ht_item Ht_item;
	 
	// Define the Hash Table Item here
	struct Ht_item {
		char* key;
	};
	 
	 
	typedef struct LinkedList LinkedList;
	 
	// Define the Linkedlist here
	struct LinkedList {
		Ht_item* item; 
		LinkedList* next;
	};
	 
	// Define the Hash Table here
	struct HashTable {
		// Contains an array of pointers
		// to items
		Ht_item** items;
		LinkedList** overflow_buckets;
		int size;
		int count;
	};

	typedef struct HashTable HashTable;

	HashTable* ht_create_table();
	void ht_insert(HashTable* table, char* key);
	void print_table(HashTable* table);
	bool ht_check(HashTable* table, char* key);
	void print_search(HashTable* table, char* key);
	
#ifdef __cplusplus
}
#endif

#endif 