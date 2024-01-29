# Generate synthetic seismograms from all SAC files in the directory
# using INSTASEIS and saves as SAC files at the output directory.
#
# It assumes an environment variable $IFILES to be a directory and the
# CMT catalog files have to be located at $IFILES/CMT/
#
# How to run from terminal:
# >>> python /path/to/getsynthetics.py /path/to/sacdir /path/to/outputdir model
# if model is not specified, it is assumed to be 'ak135f_1s'
#
# Last modified by sirawich-at-princeton.edu, 11/15/2021

import obspy
import instaseis
import sys
import os

def main():
    args = sys.argv[0:]
    print('length(args) = %d' % len(args))
    if len(args) == 4:
        print('run getsynthetics("%s", "%s", "%s")' % (args[1], args[2], args[3]))
        getsynthetics(args[1], args[2], args[3])
    elif len(args) == 3:
        print('run getsynthetics("%s", "%s")' % (args[1], args[2]))
        getsynthetics(args[1], args[2])
    else:
        print('The number of input argument must be 3 or 4.')


def getsynthetics(sacdir, outputdir, model='ak135f_1s'):
    print("successfully parsed arguments")

    # read seismograms
    st = obspy.read(sacdir + '*.sac', format='SAC')

    print("successfully read SAC files")
    print("%d files are read" % len(st))

    # read CMT catalog
    dt = st[0].stats.starttime + st[0].stats.sac['user8']
    fname = "%sCMT/%s%02d.ndk" % (os.environ.get('IFILES'), dt.ctime()[4:7].lower(), dt.year % 100)
    try:
        cat = obspy.read_events(fname)
    except FileNotFoundError:
        print("No such file or directory. Exit.\n")
        return
    
    print("successfully read the GCMT file")
    print("%d events are read" % len(cat))

    # filter CMT catalog
    lat = st[0].stats.sac['evla']
    lon = st[0].stats.sac['evlo']
    mag = st[0].stats.sac['mag']
    dep = st[0].stats.sac['evdp'] * 1000

    cat = cat.filter("magnitude >= %f" % (mag - 0.5),
                     "magnitude <= %f" % (mag + 0.5),
                     "latitude >= %f" % (lat - 1),
                     "latitude <= %f" % (lat + 1),
                     "longitude >= %f" % (lon - 1),
                     "longitude <= %f" % (lon + 1),
                     "depth >= %f" % (dep - 50000),
                     "depth <= %f" % (dep + 50000),
                     "time > %s" % str(dt - 60),
                     "time < %s" % str(dt + 60))

    print("Successfully filtered the catalog")
    print("%d events remained" % len(cat))
    print(cat)

    if len(cat) != 1:
        print("Cannot find the moment tensor for the event. Exit.\n")
        return

    tensor = cat[0].focal_mechanisms[0].moment_tensor.tensor
    source = instaseis.Source(latitude=lat, longitude=lon, depth_in_m=dep,
                              m_rr=tensor.m_rr,
                              m_tt=tensor.m_tt,
                              m_pp=tensor.m_pp,
                              m_rt=tensor.m_rt,
                              m_rp=tensor.m_rp,
                              m_tp=tensor.m_tp,
                              origin_time=dt)

    print("Successfully loaded the source")
    print(tensor)
    # read model file
    db = instaseis.open_db("syngine://%s" % model)

    print("Successfully loaded the Earth model")

    # make the output directory if it does not exist
    os.makedirs(outputdir, exist_ok=True)

    for tr in st:
        # TODO: figure out how to use receiver depth with syngine database
        #       right now instaseis does not use station depth
        rec = instaseis.Receiver(latitude=tr.stats.sac['stla'],
                                 longitude=tr.stats.sac['stlo'],
                                 depth_in_m=-tr.stats.sac['stel'],
                                 network=tr.stats.network,
                                 station=tr.stats.station)

        st_s = db.get_seismograms(source=source, receiver=rec)

        print("Getting synthetic seismograms ...")

        for ii, tr_s in enumerate(st_s):
            # copy header variables to newly created seismograms
            tr_s.stats.sac = tr.stats.sac
            # modify sac header to agree with other stats parameters
            tr_s.stats.sac['delta'] = tr_s.stats.delta
            tr_s.stats.sac['scale'] = 1
            dt_s = tr_s.stats.starttime  # start time
            dt_e = tr_s.stats.endtime  # end   time
            tr_s.stats.sac['b'] = dt_s.microsecond % 1000
            tr_s.stats.sac['e'] = tr_s.stats.sac['b'] + (dt_e - dt_s)
            tr_s.stats.sac['nzyear'] = dt_s.year
            tr_s.stats.sac['nzjday'] = dt_s.julday
            tr_s.stats.sac['nzhour'] = dt_s.hour
            tr_s.stats.sac['nzmin'] = dt_s.minute
            tr_s.stats.sac['nzsec'] = dt_s.second
            tr_s.stats.sac['nzmsec'] = dt_s.microsecond // 1000
            tr_s.stats.sac['kcmpnm'] = tr_s.stats.channel
            # remove any arrival time tags as they are no longer relevant
            tr_s.stats.sac['t0'] = -12345
            tr_s.stats.sac['kt0'] = '-12345  '
            # give a name for the sac file
            name = str(tr_s.stats.starttime)[:-8] + '_' + tr_s.stats.station[-2:]
            name = name.replace('-', '').replace(':', '')

            print(tr_s.stats)
            # write the sac file
            tr_s.write("%s%s_%d_%s_SYNTHETIC.sac" % (outputdir, name, ii, model), format="SAC")

    print("Everything is done: Exit.\n")


if __name__ == "__main__":
    main()
