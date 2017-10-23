#ifndef _GNU_SOURCE /* These are needed to work with glib */
#define _GNU_SOURCE
#endif
#include <dlfcn.h>
#include <glib.h> /* gcc flags will provide us with the write linker options */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 * This file will contain the list of blacklist functions and
 * whitelist functions that we want to search for during runtime.
 *
 * This file is generated with a search before compiling your program
 * and then referenced in the code below.
 */
#include "tracefunc.h"

static FILE* fp_trace;

static int nb_programs = 0;

/*
 * We will store the 'dlopen' handles used to probe each binary/library.
 */
static void** lookup_handles;

/*
 * This hashtable will store the addresses of the human-readable functions
 * that we will be looking for while the program is running.
 *
 * Each time a function is traced, we look it up here and spit it out if
 * we get a match table and we skip if we get a hit in the blacklist table.
 */
GHashTable *match_addrs, *black_addrs;

/*
 * The following attributes and prototypes are important:
 *
 * 'constructor' / 'destructor' is required by GCC to initiate tracing.
 *
 * 'no_instrument_function' is optional, but important, otherwise GCC
 * will trace our tracing code, which easily leads to infinite recursion
 * and will cause the program to segfault.
 */
void __attribute__((constructor, no_instrument_function)) trace_begin(void);
void __attribute__((destructor, no_instrument_function)) trace_end(void);
void __attribute__((no_instrument_function))
print(const char* direction, void* func, void* caller);
void __attribute__((no_instrument_function))
__cyg_profile_func_enter(void* func, void* caller);
void __attribute__((no_instrument_function))
__cyg_profile_func_exit(void* func, void* caller);

/*
 * Perform the translation between function name and address
 * and store it in a hashtable for later usage.
 */
void lookup_function(const char* function_name, GHashTable* table)
{
    int x;
    void* result;
    for (x = 0; x < nb_programs; x++)
    {
        result = dlsym(lookup_handles[x], function_name);
        if (result)
        {
            g_hash_table_insert(table, result, (void*)1);
            break;
        }
    }

    if (result == NULL)
    {
        printf("match function %s cannot be found in any of the libraries\n",
               function_name);
        exit(1);
    }
}

/*
 * Invoke 'dlopen' for each external library we want to query during tracing,
 * including ourselves (the main program).
 */
void load_library(const char* library, int entry)
{
    lookup_handles[entry] = dlopen(library, RTLD_NOW);
    if (lookup_handles[entry] == NULL)
    {
        printf("Could not open binary/library named: %s, because: %s\n",
               library ? library : "(main)", dlerror());
        exit(1);
    }
    printf("Opened library located at %s\n", library ? library : "(main)");
}
void __attribute__((constructor, no_instrument_function)) trace_begin(void)
{
    int match_count, black_count;
    int x = 0;

    printf("program start\n");
    /*
     * First, use 'dlopen' to open a handle to each binary/library
     * that we intend to probe while tracing the program.
     */
    nb_programs = sizeof(programs) / sizeof(programs[0]);

    printf("Will probe %d libraries, including main(), for tracing...\n",
           nb_programs + 1);

    lookup_handles = malloc((nb_programs + 1) * sizeof(void*));

    /* Open external libraries library */
    for (x = 0; x < nb_programs; x++)
        load_library(programs[x], x);

    /* Open the main program */
    load_library(NULL, nb_programs);
    nb_programs++;

    /*
     * Initialize the match and blacklist hashtables.
     */
    match_addrs = g_hash_table_new(g_direct_hash, g_direct_equal);
    black_addrs = g_hash_table_new(g_direct_hash, g_direct_equal);

    /*
     * Open the output trace results file.
     */
    fp_trace = fopen(output_trace_filename, "w");

    if (fp_trace == NULL)
    {
        perror("fopen");
        printf("Failed to open output trace file: %s\n", output_trace_filename);
    }

    match_count = sizeof(matches) / sizeof(matches[0]);
    black_count = sizeof(blacklist) / sizeof(blacklist[0]);

    /*
     * Now, go through each of the requested human-readable function names and
     * find the corresponding address of each function. Store that address
     * into the hashtables so that we can figure out whether or not to output
     * each traced function while the program is running.
     */
    printf("Looking up addresses for %d whitelisted functions...\n",
           match_count);

    for (x = 0; x < match_count; x++)
        lookup_function(matches[x], match_addrs);

    printf("Looking up addresses for %d blacklisted functions...\n",
           black_count);

    for (x = 0; x < black_count; x++)
        lookup_function(blacklist[x], black_addrs);
}
void __attribute__((destructor, no_instrument_function)) trace_end(void)
{
    if (fp_trace != NULL)
    {
        fclose(fp_trace);
    }
    printf("program end\n");
}

void __attribute__((no_instrument_function))
print(const char* direction, void* func, void* caller)
{
    Dl_info dl1, dl2;

    if (fp_trace == NULL)
        return;

    /*
     * Did the currently traced function 'hit' in the match table?
     */
    if (!g_hash_table_lookup(match_addrs, func))
        return;

    /*
     * Match hit. Check the blacklist.
     * Did we hit in the blacklist? Then ignore this one.
     */
    if (g_hash_table_lookup(black_addrs, func))
        return;

    /*
     * Now that we know the addresses in question are found,
     * we need to print out human-readable results by converting
     * the addresses of both the caller and callee to function names.
     */
    dladdr(func, &dl1);
    dladdr(caller, &dl2);

    /*
     * Sometimes this happens, no idea why.
     */
    if (dl1.dli_sname == NULL)
        return;

    /* Finished. */
    fprintf(fp_trace, "time [%ld] addr (%p): %s call from (%s) => to (%s) \n",
            time(NULL), func, direction,
            dl2.dli_sname ? dl2.dli_sname : "unknown",
            dl1.dli_sname ? dl1.dli_sname : "unknown");
    fflush(fp_trace);
}

/*
 * This functions are required to be defined by GCC.
 * Each traced function results in GCC invoking these functions,
 * from which we do our more sophisticated tracing.
 */
void __attribute__((no_instrument_function))
__cyg_profile_func_enter(void* func, void* caller)
{
    print("enter", func, caller);
}

void __attribute__((no_instrument_function))
__cyg_profile_func_exit(void* func, void* caller)
{
    print("exit", func, caller);
}
