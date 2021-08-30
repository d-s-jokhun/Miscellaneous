

#%%

from os import listdir, rename
from os.path import isfile, join
import multiprocessing as mp
import shutil

#%%


src_dir = "/MBELab/jokhun/Pro 1/Cellular biosensor/ML temp/raw stacks"
dst_dir = "/MBELab/jokhun/Pro 1/Cellular biosensor/ML temp/projected"

# id_str = 'Med_Sam2__w1CSU DIC_'
# replacement_str = 'Med_Fea1_Sam2_DIC_'
# filenames = [f for f in listdir(src_dir) if isfile(join(src_dir, f)) and id_str in f]

id_str = '_w1CSU DIC_'
replacement_str = '_DIC_'
filenames = [f for f in listdir(src_dir) if isfile(join(src_dir, f)) and id_str in f and 'Med_Fea' in f]

src_paths = [join(src_dir,f) for f in filenames]
dst_paths = [join(dst_dir,f.replace(id_str,replacement_str)[:-4]+".TIF") for f in filenames]


#%%

if __name__ == '__main__':
    with mp.Pool() as p:
        p.starmap (shutil.copyfile,zip(src_paths,dst_paths))



#%%




#%%




