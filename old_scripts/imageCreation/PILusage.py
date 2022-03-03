from PIL import Image
import numpy as np

def random_img(output, width, height):

    array = np.randint(0,255, (height,width,3))  

    array = np.array(array, dtype=np.uint8)
    img = Image.fromarray(array)
    img.save(output)


random_img('random.png', 100, 50)