'''
Inference

This script infers the aesthetic value of all images. If several weights are 
available the most recent one will be chosen 


authors : Juliette Langlois <juliette.a.langlois@gmail.com>
          Nicolas Mouquet <nicolas.mouquet@cnrs.fr>
          Cedric Braga <cedric.braga@hotmail.fr>
          Valentine Fleuré <valentine.fleure@gmail.com>

date : 2021/07/30
'''

#Misc librairies
import os
import glob
import multiprocessing
import numpy as np
import pandas as pd
from tqdm import tqdm
from PIL import Image
from pathlib import Path
import torch
from torch import nn
from torch.utils.data import Dataset
from torchvision import models, transforms
import shutil

multiprocessing.cpu_count()

class EsthDataset(torch.utils.data.Dataset):
    def __init__(self, paths, transforms=None):
        self.paths = paths
        self.transform = transforms
        
    def __len__(self):
        return len(self.paths)
    
    def __getitem__(self, index):
        image = Image.open(self.paths[index]).convert('RGB')
        
        if image.size[0] != 500 or image.size[1] != 500 :
            print(self.paths[index], ' is the wrong size : ', image.size, ". It has been resized to (500,500) but the accuracy of the score cannot be guaranteed")
            image = image.resize((500,500))
        
        if self.transform is not None:
            image = self.transform(image)
            
        return image, torch.tensor(index, dtype=torch.float32)

# hyperparameters config for the model architecture and data loader

# Use GPU or CPU
use_cuda = torch.cuda.is_available()
device = torch.device("cpu")

run = '1'
#Path for the weights of the saved model we want to use
# Latest weights
list_of_files = glob.glob(('../data/raw_data/aesthetic_deep_model_Langlois_2022/03_train_all/*.pth'))
print(list_of_files)
weights_path =  max(list_of_files, key=os.path.getctime)

#Path to the input data and output file
in_path = Path('../data/raw_data/aesthetic_deep_model_Langlois_2022/png/')
out_path = Path('../outputs/2g_aesthetic_inference_new_sp.csv')

def get_model_for_eval(path_to_weights):
    """Gets the broadcasted model."""
    model = models.resnet50()
    
    model.fc = nn.Sequential(
    nn.Dropout(0.5),
    nn.Linear(2048, 1)
    )
    
    model.load_state_dict(torch.load(path_to_weights))
    model.eval()
    
    return model

def predict_batch(paths, path_to_weights, output_path):
    transform = transforms.Compose([
    transforms.Resize((224,224)),
    transforms.ToTensor(),
    transforms.Normalize([0.865, 0.877, 0.881], [0.273, 0.252, 0.250])
    ])
    
    model = get_model_for_eval(path_to_weights)
    model.to(device)
    
    images = EsthDataset(paths, transforms=transform)
    loader = torch.utils.data.DataLoader(images, batch_size=4, shuffle=False)
    
    im_names = torch.tensor([], dtype=torch.float, device=device)
    y_pred = torch.tensor([], device=device)
    
    with torch.no_grad():
        for data in tqdm(loader):
            inputs = [i.to(device) for i in data[:-1]]
            names = data[-1].to(device)
            outputs = model(*inputs)
            im_names = torch.cat((im_names, names), 0)
            y_pred = torch.cat((y_pred, outputs), 0)
    y_pred = y_pred.cpu().numpy().flatten()
    im_names = im_names.cpu().numpy().flatten().astype(int)
    df = pd.DataFrame({"image_name":[paths[i].split('/')[-1] for i in im_names],"predicted_score":y_pred})
    df.to_csv(output_path, sep=",", index=False)

files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(in_path) for f in filenames if os.path.splitext(f)[1] in ['.png', '.jpg','.jpeg']]
predict_batch(files, weights_path, out_path)
