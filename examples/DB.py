#!/usr/bin/env python


"""
Some simple functions for use with user DBs in the DES DB.
"""

import suchyta_utils as es


if __name__ == "__main__":
  
    # Find all your tables in the DB
    # You can replace None with the name of another user to check their DB tables, e.g. SUCHYTA1
    tables = es.db.GetTableNames(user=None)
    print 'found %i tables' %(len(tables))


    # Search for any tables that have the string 'DEBUG' in them
    dbg = es.db.SearchTables(tables, 'DBG')
    print dbg


    '''
    # Delete these tables.
    # (Don't worry, others don't have permissions to delete your tables.)
    es.db.Drop(dbg)
    '''


    #Check your DB quota
    print es.db.GetQuota()


    # Print a summary of your DB usage
    es.db.PrintUsage(user='JOEZUNTZ')
