"""
The :mod:`suchyta_utils.db` submodule is for working with the DES Oracle database, particularly the user DB space.
Mostly, it's functions for calling Oracle commands whose syntax I tend to forget,
that I don't want to have to google every time I use them.
(Like most physicists, I'm far from an Oracle expert.)


.. note::
    For :mod:`suchyta_utils.db` to work, you need to have your .netrc file setup for the DESDM DB.
    (see the `desdb github repository <https://github.com/esheldon/desdb#access-to-servers>`_).

"""

import _dbfunctions
import desdb as _desdb
import numpy as _np



def _get_table_names():
    cur = _dbfunctions.get_cursor()
    cur.execute("select owner, table_name from dba_tables order by table_name")
    return _np.array( cur.fetchall() )[:,:]


def GetTableNames(user=None):
    """
    Get all of a user's tables from the DB.

    Parameters
    ----------
    user (str)
        The username whose tables you want to return. None means use your personal username (read from the .netrc file).

    Returns
    -------
    tables (list)
        List of the table's owned by the specified user

    """
    if user is None:
        user, pwd = _dbfunctions.retrieve_login(_dbfunctions.db_specs.db_host)
        user = user.upper()

    names = _get_table_names()
    tables = []
    for name in names:
        if name[0].find(user)!=-1:
            tables.append(name[1])
    return tables


def Drop(tables):
    """
    Delete tables from the DB. (You must own them.)

    Parameters
    ----------
    tables (list)
        A list, where each entry is a (str) table name

    Returns
    -------
    None

    """
    cur = _desdb.connect()
    for table in tables:
        print table
        cur.quick("DROP TABLE %s PURGE" %table)
    cur.commit()


def GetQuota(user=None):
    """
    Get a user's DB quota (in GB).

    Parameters
    ----------
    user (str)
        The username whose quota you want to return. None means use your personal username (read from the .netrc file).

    Returns
    -------
    quota (float)
        The user's quota in GB.

    """
    if user is None:
        user, pwd = _dbfunctions.retrieve_login(_dbfunctions.db_specs.db_host)
        user = user.upper()

    cur = _desdb.connect()
    q = "SELECT USERNAME, TABLESPACE_NAME, MAX_BYTES from DBA_TS_QUOTAS WHERE USERNAME='%s'" %(user)
    all = cur.quick(q, array=True)
    return all['max_bytes'][0] / _np.power(1024., 3)


def _add(tables, user='SUCHYTA1'):
    cur = _desdb.connect()
    count = 0
    for n in tables:
        q = """SELECT SUM(bytes), SUM(bytes) B from dba_extents where owner='%s' and segment_name='%s'""" %(user, n)

        '''
        all = cur.quick(q, array=True)
        count += all['b'][0]
        '''
        try:
            all = cur.quick(q, array=True)
            count += all['b'][0]
            print n, '  ', all['b'][0] / _np.power(1024.0, 3), 'GB'
        except:
            pass

    print count/_np.power(1024.0, 3), 'GB'


def TotalUsage(user=None):
    """
    Find total DB usage for the user.

    Parameters
    ----------
    user (str)
        The username whose usage you want to return. None means use your personal username (read from the .netrc file).

    Returns
    -------
    usage (float)
        Usage in GB

    """

    if user is None:
        user, pwd = _dbfunctions.retrieve_login(_dbfunctions.db_specs.db_host)
        user = user.upper()

    cur = _desdb.connect()
    return cur.quick("""SELECT SUM(bytes) as b from dba_extents where owner='%s'"""%(user), array=True)['b'][0] / _np.power(1024.0,3)


def PrintUsage(user=None):
    """
    Print a summary of a user's DB usage.
    Each line prints a table, and its GB usage.
    The last line prints the total usage

    Parameters
    ----------
    user (str)
        The username whose quota you want to return. None means use your personal username (read from the .netrc file).

    Returns
    -------
    None

    """
    if user is None:
        user, pwd = _dbfunctions.retrieve_login(_dbfunctions.db_specs.db_host)
        user = user.upper()

    tables = GetTableNames(user=user)
    _add(tables, user=user)


def SearchTables(tables, key):
    """
    Find which tables in a list contain a string.

    Parameters
    ----------
    tables (list)
        A list, where each entry is a (str) table name
    key (str)
        String to search the list for

    Returns
    -------
    matches (list)
        A (str) list of the tables identified as matches

    """
    ts = []
    for name in tables:
        if name.find(key)!=-1:
            ts.append(name)
    return ts


def ColumnDescribe(table, user=None):
    """
    Get the column names, as well as data type and precsion information about a DB table

    Parameters
    ----------
    table (str)
        Name of the table
    user (str)
        Name of the user who owns the table. None does not specify a user. (A user is not needed if the table name is unique.)

    Returns
    -------
    arr (structured array)
        Array with columns [column_name, data_type, data_precision, data_scale]

    """
    cur = _desdb.connect()
    if user is not None:
        arr = cur.quick("SELECT column_name, data_type, data_precision, data_scale from all_tab_cols where table_name='%s' and owner='%s'"%(table.upper(),user.upper()), array=True)
    else:
        arr = cur.quick("SELECT column_name, data_type, data_precision, data_scale from all_tab_cols where table_name='%s'"%(table.upper()), array=True)
    return arr


def IndexDescribe(table, user=None):
    """
    Get information about the indexing of a table. 

    Parameters
    ----------
    table (str)
        Name of the table
    user (str)
        Name of the user who owns the table. None does not specify a user. (A user is not needed if the table name is unique.)

    Returns
    -------
    arr (structured array)
        Array with several columns [index_owner, index_name, table_owner, table_name, column_name, column_position, column_length, char_length, descend] 

    """
    cur = _desdb.connect()
    if user is not None:
        arr = cur.quick("select * from dba_ind_columns where table_name='%s' and table_owner='%s'"%(table.upper(),user.upper()), array=True)
    else:
        arr = cur.quick("select * from dba_ind_columns where table_name='%s'"%(table.upper()), array=True)
    return arr


"""
    For indexing Balrog tables. Don't use this if you don't know what you're doing.
"""
def IndexBalrog(db, dname, tab, what, name):
    cur = _desdb.connect()

    bands = ['det', 'g', 'r', 'i', 'z', 'y']
    for band in bands:
        q = "CREATE INDEX i_%s_%s%s_%s ON balrog_%s_%s_%s %s" %(dname,tab[0],band[0],name, db,tab,band,what)
        print q
        arr = cur.quick(q, array=True)


'''
def print_table_cols(table):
    cur = _dbfunctions.get_cursor()
    cur.execute("select * from %s" %table)
    tab = _np.array( cur.description )[:,0] 
    for t in tab:
        print t


def get_table_cols(table):
    cur = _dbfunctions.get_cursor()
    cur.execute("select * from %s" %table)
    tab = _np.array( cur.description )[:,0] 
    return tab


def test_entry(table):
    cur = _dbfunctions.get_cursor()
    cur.execute("select * from %s" %table)
    return cur.fetchone()
'''

'''
def GetHealPixRectangles(nside, index, nest=False):
    vec_corners = hp.boundaries(nside, index, nest=nest)
    vec_corners = _np.transpose(vec_corners, (0,2,1))
    vec_corners = _np.reshape(vec_corners, (vec_corners.shape[0]*vec_corners.shape[1], vec_corners.shape[2]))
   
    theta_corners, phi_corners = hp.vec2ang(vec_corners)
    theta_corners = _np.reshape(theta_corners, (theta_corners.shape[0]/4, 4))
    phi_corners = _np.reshape(phi_corners, (phi_corners.shape[0]/4, 4))

    ra_corners = _np.degrees(phi_corners)
    dec_corners = 90.0 - _np.degrees(theta_corners)

    ramin = _np.amin(ra_corners, axis=-1)
    ramax = _np.amax(ra_corners, axis=-1)
    decmin = __np.amin(dec_corners, axis=-1)
    decmax = _np.amax(dec_corners, axis=-1)

    return ramin, ramax, decmin, decmax

'''
