#!/usr/bin/env python


"""
Some simple functions for use with user DBs in the DES DB.
"""

import suchyta_utils as es


if __name__ == "__main__":
  
    # Find all your tables in the DB
    # You can replace None with the name of another user to check their DB tables, e.g. SUCHYTA1
    tables = es.db.get_my_tables(user=None)


    # Search for any tables that have the string 'DEBUG' in them
    dbg = es.db.search_tables(tables, 'DEBUG')
    print dbg


    # Delete these tables.
    # (Don't worry, others don't have permissions to delete your tables.)
    es.db.drop(dbg)


    #Check your DB quota
    print es.db.get_quota()


    # Print a summary of your DB usage
    es.db.check_usage()
