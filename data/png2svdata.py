import numpy as np
import cv2
import os

src_folder_name = 'dataset/ori/mito'
dst_folder_name = 'dataset/processed/mito'

os.makedirs("dataset/processed/mito", exist_ok=True)

dim_z = len(os.listdir(src_folder_name))

for file_index in range(0, dim_z):
    array = cv2.imread(f"{src_folder_name}/{file_index:04d}.png")[:, :, 0]
    array = np.array(array, dtype=np.int)

    with open(f"{dst_folder_name}/{file_index:04d}.svdata", "wb") as f_out:
        sb = (int(array.shape[0])).to_bytes(4, byteorder='little')
        f_out.write(sb)
        sb = (int(array.shape[1])).to_bytes(4, byteorder='little')
        f_out.write(sb)
        sb = (4).to_bytes(4, byteorder='little')
        f_out.write(sb)

        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                sb = int(array[i, j])
                sb = sb.to_bytes(4, byteorder='little')
                f_out.write(sb)
