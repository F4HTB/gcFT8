// Open-addressing callsign set used by the ADIF already-worked filter.

#include "hash_table.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HASH_INITIAL_CAPACITY 1024u
#define HASH_MAX_LOAD_NUMERATOR 7u
#define HASH_MAX_LOAD_DENOMINATOR 10u

static uint32_t hash_function(const char* str)
{
    uint32_t hash = 2166136261u;

    while (*str)
    {
        hash ^= (unsigned char)*str;
        hash *= 16777619u;
        ++str;
    }

    return hash;
}

static bool normalize_key(const char* key, char normalized[HASH_TABLE_KEY_SIZE])
{
    size_t start = 0;
    size_t end;
    size_t out_len = 0;

    if ((key == NULL) || (normalized == NULL))
        return false;

    while ((key[start] != '\0') && isspace((unsigned char)key[start]))
        ++start;

    end = strlen(key);
    while ((end > start) && isspace((unsigned char)key[end - 1]))
        --end;

    if (end == start)
        return false;

    for (size_t idx = start; idx < end; ++idx)
    {
        if ((out_len + 1) >= HASH_TABLE_KEY_SIZE)
            return false;

        normalized[out_len++] = (char)toupper((unsigned char)key[idx]);
    }

    normalized[out_len] = '\0';
    return true;
}

static size_t next_power_of_two(size_t value)
{
    size_t result = 1;

    while (result < value)
        result <<= 1;

    return result;
}

static bool ht_find_slot(const HashTable* table, const char* key, uint32_t hash, size_t* slot, bool* found)
{
    size_t mask;
    size_t idx;

    if ((table == NULL) || (table->items == NULL) || (table->size == 0) || (slot == NULL) || (found == NULL))
        return false;

    mask = table->size - 1u;
    idx = (size_t)hash & mask;

    for (size_t probe = 0; probe < table->size; ++probe)
    {
        Ht_item* item = &table->items[idx];

        if (!item->used)
        {
            *slot = idx;
            *found = false;
            return true;
        }

        if ((item->hash == hash) && (strcmp(item->key, key) == 0))
        {
            *slot = idx;
            *found = true;
            return true;
        }

        idx = (idx + 1u) & mask;
    }

    return false;
}

static bool ht_init_items(HashTable* table, size_t capacity)
{
    capacity = next_power_of_two(capacity);
    table->items = (Ht_item*)calloc(capacity, sizeof(table->items[0]));
    if (table->items == NULL)
        return false;

    table->size = capacity;
    table->count = 0;
    return true;
}

static bool ht_insert_normalized(HashTable* table, const char* key, uint32_t hash, bool warn_duplicate)
{
    size_t slot = 0;
    bool found = false;

    if (!ht_find_slot(table, key, hash, &slot, &found))
        return false;

    if (found)
    {
        if (warn_duplicate)
            printf("Attention, duplicate entry in QSO filter table: %s\n", key);
        return true;
    }

    table->items[slot].used = true;
    table->items[slot].hash = hash;
    strcpy(table->items[slot].key, key);
    ++table->count;
    return true;
}

static bool ht_grow(HashTable* table)
{
    Ht_item* old_items;
    size_t old_size;
    size_t old_count;

    if ((table == NULL) || (table->size > (SIZE_MAX / 2u)))
        return false;

    old_items = table->items;
    old_size = table->size;
    old_count = table->count;

    table->items = NULL;
    table->size = 0;
    table->count = 0;

    if (!ht_init_items(table, old_size * 2u))
    {
        table->items = old_items;
        table->size = old_size;
        table->count = old_count;
        return false;
    }

    for (size_t idx = 0; idx < old_size; ++idx)
    {
        if (old_items[idx].used)
            (void)ht_insert_normalized(table, old_items[idx].key, old_items[idx].hash, false);
    }

    free(old_items);
    return true;
}

HashTable* ht_create_table(void)
{
    HashTable* table = (HashTable*)malloc(sizeof(*table));
    if (table == NULL)
        return NULL;

    table->items = NULL;
    table->size = 0;
    table->count = 0;

    if (!ht_init_items(table, HASH_INITIAL_CAPACITY))
    {
        free(table);
        return NULL;
    }

    return table;
}

void free_table(HashTable* table)
{
    if (table == NULL)
        return;

    free(table->items);
    free(table);
}

void ht_insert(HashTable* table, const char* key)
{
    char normalized[HASH_TABLE_KEY_SIZE];
    uint32_t hash;

    if ((table == NULL) || !normalize_key(key, normalized))
        return;

    if (((table->count + 1u) * HASH_MAX_LOAD_DENOMINATOR) >= (table->size * HASH_MAX_LOAD_NUMERATOR))
    {
        if (!ht_grow(table))
        {
            printf("Insert Error: Hash Table is full\n");
            return;
        }
    }

    hash = hash_function(normalized);
    if (!ht_insert_normalized(table, normalized, hash, true))
        printf("Insert Error: Hash Table is full\n");
}

bool ht_check(const HashTable* table, const char* key)
{
    char normalized[HASH_TABLE_KEY_SIZE];
    uint32_t hash;
    size_t slot = 0;
    bool found = false;

    if ((table == NULL) || !normalize_key(key, normalized))
        return false;

    hash = hash_function(normalized);
    if (!ht_find_slot(table, normalized, hash, &slot, &found))
        return false;

    return found;
}

static const char* ht_search(const HashTable* table, const char* key)
{
    char normalized[HASH_TABLE_KEY_SIZE];
    uint32_t hash;
    size_t slot = 0;
    bool found = false;

    if ((table == NULL) || !normalize_key(key, normalized))
        return NULL;

    hash = hash_function(normalized);
    if (!ht_find_slot(table, normalized, hash, &slot, &found) || !found)
        return NULL;

    return table->items[slot].key;
}

void print_search(const HashTable* table, const char* key)
{
    const char* val = ht_search(table, key);
    if (val == NULL)
    {
        printf("%s does not exist\n", key);
        return;
    }

    printf("Key:%s Value:%s\n", key, val);
}

void print_table(const HashTable* table)
{
    if (table == NULL)
        return;

    printf("\n-------------------\n");
    for (size_t idx = 0; idx < table->size; ++idx)
    {
        if (table->items[idx].used)
            printf("Index:%zu, Key:%s\n", idx, table->items[idx].key);
    }
    printf("-------------------\n");
}
