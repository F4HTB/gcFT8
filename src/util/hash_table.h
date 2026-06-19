#ifndef _HASH_H_
#define _HASH_H_

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define HASH_TABLE_KEY_SIZE 20

#ifdef __cplusplus
extern "C"
{
#endif

	typedef struct Ht_item Ht_item;

	struct Ht_item {
		char key[HASH_TABLE_KEY_SIZE];
		uint32_t hash;
		bool used;
	};

	struct HashTable {
		Ht_item* items;
		size_t size;
		size_t count;
	};

	typedef struct HashTable HashTable;

	HashTable* ht_create_table(void);
	void free_table(HashTable* table);
	void ht_insert(HashTable* table, const char* key);
	bool ht_check(const HashTable* table, const char* key);
	void print_table(const HashTable* table);
	void print_search(const HashTable* table, const char* key);

#ifdef __cplusplus
}
#endif

#endif
