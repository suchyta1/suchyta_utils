"""
Test Stuff

"""

import _dbfunctions
import desdb as _desdb
import numpy as _np



def _get_table_names():
    cur = _dbfunctions.get_cursor()
    cur.execute("select owner, table_name from dba_tables order by table_name")
    return _np.array( cur.fetchall() )[:,:]


def get_my_tables(user=None):
    """
    Return all of a user's  tables from the DB.
    Default user=None returns your own tables.

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

def drop(tables):
    """
    Delete the tables in the list. You must own them.

    """
    cur = _desdb.connect()
    for table in tables:
        print table
        cur.quick("DROP TABLE %s PURGE" %table)
    cur.commit()


def get_quota(user=None):
    """
    Return user's DB quota.
    Default user=None returns your own quota.

    """
    if user is None:
        user, pwd = _dbfunctions.retrieve_login(_dbfunctions.db_specs.db_host)
        user = user.upper()

    cur = _desdb.connect()
    q = "SELECT USERNAME, TABLESPACE_NAME, MAX_BYTES from DBA_TS_QUOTAS WHERE USERNAME='%s'" %(user)
    all = cur.quick(q, array=True)
    return all['max_bytes'][0] / _np.power(1024., 3), 'GB'


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

def check_usage(user=None):
    """
    Print a summary of a user's DB usage
    Default user=None prints your own usage

    """
    if user is None:
        user, pwd = _dbfunctions.retrieve_login(_dbfunctions.db_specs.db_host)
        user = user.upper()

    tables = get_my_tables(user=user)
    _add(tables, user=user)


def search_tables(tables, key):
    """
    Search a list of tables for a string, and return any containing that string

    Arguments:
    tables: a list of table names
    key: string to search for

    Returns:
    a list of the tables identified

    """
    ts = []
    for name in tables:
        if name.find(key)!=-1:
            ts.append(name)
    return ts


'''
def Cat2New(t1, t2, outtab, inc):
    cur = _desdb.connect()

    q = """BEGIN
            EXECUTE IMMEDIATE 'DROP TABLE %s'; \
            EXCEPTION
                WHEN OTHERS THEN \
                    IF SQLCODE != -942 THEN \
                        RAISE; \
                    END IF; \
            END;""" %(outtab)
    cur.quick(q, array=True)

    q = "CREATE TABLE %s AS SELECT * from %s" %(outtab, t1)
    cur.quick(q)
    info = cur.quick("select column_name, data_type, data_length from all_tab_columns where upper(table_name) = '%s'"%(t2), array=True)
    outnames_i = ', '.join(info['column_name'])
    outnames_s = outnames_i.replace('BALROG_INDEX', 'BALROG_INDEX + %f'%inc)
    q = """INSERT INTO %s (%s) SELECT %s FROM %s""" %(outtab, outnames_i, outnames_s, t2)
    cur.quick(q) 
    cur.quick('COMMIT WORK')


def Combine(r1, r2, outr, bands=['G','R','I','Z','Y'], types=['TRUTH','NOSIM','SIM','DES']):
    cur = _desdb.connect()
    all = cur.quick("select count(*) as count from %s_TRUTH_G" %(r1), array=True)
    print all
    inc = all['count'][0]

    for band in bands:
        for type in types:
            outname = '%s_%s_%s' %(outr, type, band)
            rr1 = '%s_%s_%s' %(r1, type, band)
            rr2 = '%s_%s_%s' %(r2, type, band)
            print rr1, rr2, '-->', outname
            Cat2New(rr1, rr2, outname, inc)
'''


def IndexBalrog(db, dname, tab, what, name):
    """
    For indexing Balrog tables. Don't use this if you don't know what you're doing.
    """
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



def Cat2Existing(existing, new):
    cur = _desdb.connect()
    info = cur.quick("select count(*) as num from %s" %(existing), array=True)
    num = int(info['num'])
    
    info = cur.quick("select column_name, data_type, data_length from all_tab_columns where upper(table_name) = '%s'"%(new), array=True)
    if type!='DES':
        cut = (info['column_name']=='balrog_index')
        info['column_name'][cut] = 'balrog_index+%i' %(num)

    outnames = ', '.join(info['column_name'])
    innames = _np.core.defchararray.add(_np.array( ['input.']*len(info)), info['column_name'])
    innames = ', '.join(innames)
    q = """INSERT INTO %s (%s) 
                SELECT %s
                FROM %s input""" %(existing, outnames, innames, new)
    print q
    cur.quick(q, array=True)



In [155]: def plotband(band, num):
    bins = _np.arange(0.1, 20, 0.2); cs = 0; fw = 20; s2n = 0; deb = 11
    sim = cur.quick("select sim.flux_auto, sim.fluxerr_auto, sim.fwhm_image, sim.class_star from balrog_debug%i_sim_%s sim where sim.class_star > %f and sim.fwhm_image < %f" %(deb, band,cs,fw), array=True); print len(sim); sn = sim['flux_auto']/sim['fluxerr_auto']; cut = (sn > s2n); sim = sim[cut];
    des = cur.quick("select flux_auto, fluxerr_auto, fwhm_image from balrog_debug%i_des_%s where class_star > %f and fwhm_image < %f" %(deb, band,cs,fw), array=True); print len(des); sn = des['flux_auto']/des['fluxerr_auto']; cut = (sn > s2n); des = des[cut];
    ax = fig.add_subplot(2,2, num)
    ax.hist(sim['fwhm_image'], bins=bins, normed=True, label='sim', color='red', alpha=0.5)
    ax.hist(des['fwhm_image'], bins=bins, normed=True, label='des', color='blue', alpha=0.5)
    ax.set_xlabel('fwhm_image')
    ax.legend(loc='best')
    ax.set_title(band)
    print _np.average(sim['fwhm_image']), _np.average(des['fwhm_image'])
   .....: 

In [156]: fig = plt.figure(1, figsize=(12,7));  plotband('det', 1); plotband('r', 2); plotband('i', 3); plotband('z', 4); plt.tight_layout(); plt.show()
'''
