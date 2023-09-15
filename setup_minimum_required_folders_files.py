import os
from pathlib import Path
import urllib.request
from urllib.parse import urlparse
import ssl


tmp_file_urls = ["https://www.dropbox.com/s/bctrdr7de28m99o/B4C_Mo_Al_1-11000A.txt?dl=1", 
                 "https://www.dropbox.com/s/f51fep2nu1vr7ai/euv_sim_300_171A.fits?dl=1", 
                 "https://www.dropbox.com/scl/fi/24ckdxvwpyolfqp4120um/radiance_maps_044.fits?rlkey=pqna5na5lhj7dp3l1fe9mdq8k&dl=1",
                 "https://www.dropbox.com/scl/fi/rbgyhy7lmtftge9m8gaoq/em_map_200.sav?rlkey=4fg96y3edldbja6qy9jhp9bt0&dl=1",
                 "https://www.dropbox.com/s/z86h2h7l8pgbhnl/filter_entrance_transmission.csv?dl=1",
                 "https://www.dropbox.com/s/muclu9kncl7xyff/filter_focal_plane_transmission.csv?dl=1",
                 "https://www.dropbox.com/s/tvqt8ybipa1bm8z/aia_V9_fullemiss.nc?dl=1", 
                 "https://www.dropbox.com/s/k1qv0lujqmpvszn/si_qe_henke.csv?dl=1",
                 "https://www.dropbox.com/scl/fi/rbe7vm3sha9mbloek1iio/suncet_metadata_definition_v1.0.0.csv?rlkey=mswa2lvdrvbb9o1rer1z60p2x&dl=1"] # dl=1 is important

def run():
    if os.getenv('suncet_data') == None:
        os.environ['suncet_data'] = './'
    
    emissivity_path = Path(os.getenv('suncet_data') + '/ancillary/emissivity')
    emissivity_path.mkdir(parents=True, exist_ok=True)
    rendered_path = Path(os.getenv('suncet_data') + '/mhd/bright_fast/rendered_euv_maps')
    rendered_path.mkdir(parents=True, exist_ok=True)
    em_path = Path(os.getenv('suncet_data') + '/mhd/bright_fast/em_maps')
    em_path.mkdir(parents=True, exist_ok=True)
    reflectivity_path = Path(os.getenv('suncet_data') + '/mirror_reflectivity')
    reflectivity_path.mkdir(parents=True, exist_ok=True)
    filter_path = Path(os.getenv('suncet_data') + '/filter_transmission')
    filter_path.mkdir(parents=True, exist_ok=True)
    qe_path = Path(os.getenv('suncet_data') + '/quantum_efficiency')
    qe_path.mkdir(parents=True, exist_ok=True)
    synthetic_path = Path(os.getenv('suncet_data') + '/synthetic/level0_raw/fits')
    synthetic_path.mkdir(parents=True, exist_ok=True)
    synthetic_path = Path(os.getenv('suncet_data') + '/synthetic/level0_raw/binary')
    synthetic_path.mkdir(parents=True, exist_ok=True)
    metadata_path = Path(os.getenv('suncet_data') + '/metadata')
    metadata_path.mkdir(parents=True, exist_ok=True)

    ssl._create_default_https_context = ssl._create_unverified_context

    for url in tmp_file_urls:
        thisurl = urllib.request.urlopen(url)  
        data = thisurl.read()
        thisurl.close()
        filename = get_filename_from_url(url)
        if filename.endswith('fullemiss.nc'):
            filename = emissivity_path / filename
        elif filename.startswith('radiance_maps_'):
            filename = rendered_path / filename
        elif filename.startswith('em_map_'):
            filename = em_path / filename
        elif filename.endswith('A.txt'):
            filename = reflectivity_path / filename
        elif filename.startswith('filter_'):
            filename = filter_path / filename
        elif filename.startswith('si_qe_'):
            filename = qe_path / filename
        elif filename.startswith('suncet_metadata_'):
            filename = metadata_path / filename
        with open(filename, "wb") as f:
            print('downloading file: {}'.format(filename))
            f.write(data)

def get_filename_from_url(url):
    parsed_url = urlparse(url)
    return os.path.basename(parsed_url.path)


if __name__ == "__main__":
    run()
