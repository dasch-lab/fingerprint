# Standard imports:
import numpy as np
import torch
from torch.utils.tensorboard import SummaryWriter
from torch.utils.data import random_split
#from torch_geometric.data import DataLoader
from torch_geometric.loader import DataLoader
from torch_geometric.transforms import Compose
from pathlib import Path

# Custom data loader and model:
from data import ProteinPairsSurfaces, CenterPairAtoms
from data import RandomRotationPairAtoms, NormalizeChemFeatures, iface_valid_filter
from model import dMaSIF
#from model_copy import dMaSIF_old
from data_iteration import iterate, iterate_surface_precompute
from helper import *
from Arguments import parser
import torch.nn as nn
from pytorch_earlystopping import EarlyStopping

# Parse the arguments, prepare the TensorBoard writer:
args = parser.parse_args()
writer = SummaryWriter("runs/{}".format(args.experiment_name))
model_path = "models/" + args.experiment_name

print(f"Experiment name: {args.experiment_name} | Antibody: {args.antibody} | Flexibility: {args.flexibility}")
#print(f"Recurrent: {args.recurrent} | Weighted: {args.weighted} | Embedding dimension: {args.emb_dims}")

if not Path("models/").exists():
    Path("models/").mkdir(exist_ok=False)

# Ensure reproducibility:
torch.backends.cudnn.deterministic = True
torch.manual_seed(args.seed)
torch.cuda.manual_seed_all(args.seed)
np.random.seed(args.seed)

print(f'Flexibility: {args.flexibility} | Antibody: {args.antibody}')

# Create the model, with a warm restart if applicable:
if args.flexibility:
    args.in_channels = args.in_channels + 1
net = dMaSIF(args)
#else:
#    net = dMaSIF_old(args)
net.to(args.device)

# We load the train and test datasets.
# Random transforms, to ensure that no network/baseline overfits on pose parameters:
#args.random_rotation = True
transformations = (
    Compose([NormalizeChemFeatures(), CenterPairAtoms(), RandomRotationPairAtoms()])
    if args.random_rotation
    else Compose([NormalizeChemFeatures()])
)

# PyTorch geometric expects an explicit list of "batched variables":
batch_vars = ["xyz_p1", "xyz_p2", "atom_coords_p1", "atom_coords_p2"]
# Load the train dataset:

if args.antibody and args.mix:
    folder_data = "surface_data_ab"
elif args.antibody:
    if args.seed == 42:
        #folder_data = "surface_data_ab_pp"
        folder_data = "surface_data_42"
        print('Seed: 42')
    elif args.seed == 123:
        folder_data = "surface_data_123"
        print('Seed: 123')
    elif args.seed == 98765:
        folder_data = "surface_data_98765"
        print('Seed: 98765')
    else:
        folder_data = "surface_data_2024"
        print('Seed: 2024')
else:
    folder_data = "surface_data"

train_dataset = ProteinPairsSurfaces(
            folder_data, ppi=args.search, train=True, transform=transformations, flexibility = args.flexibility, antibody = args.antibody, mix = args.mix
        )

train_dataset = [data for data in train_dataset if iface_valid_filter(data)]
train_loader = DataLoader(
    train_dataset, batch_size=1, follow_batch=batch_vars, shuffle=True
)
print("Preprocessing training dataset")
train_dataset = iterate_surface_precompute(train_loader, net, args)

# Train/Validation split:
train_nsamples = len(train_dataset)
val_nsamples = int(train_nsamples * args.validation_fraction)
train_nsamples = train_nsamples - val_nsamples
train_dataset, val_dataset = random_split(
    train_dataset, [train_nsamples, val_nsamples]
)
print(f'Training data: {train_nsamples}')
print(f'Validation data: {val_nsamples}')

test_dataset = ProteinPairsSurfaces(
            folder_data, ppi=args.search, train=False, transform=transformations, flexibility = args.flexibility, mix = args.mix
        )


# Load the test dataset:
test_dataset = [data for data in test_dataset if iface_valid_filter(data)]
print(f'Test data: {len(test_dataset)}')
test_loader = DataLoader(
    test_dataset, batch_size=1, follow_batch=batch_vars, shuffle=True
)
print("Preprocessing testing dataset")
test_dataset = iterate_surface_precompute(test_loader, net, args)


# PyTorch_geometric data loaders:
train_loader = DataLoader(
    train_dataset, batch_size=1, follow_batch=batch_vars, shuffle=True
)
val_loader = DataLoader(val_dataset, batch_size=1, follow_batch=batch_vars)
test_loader = DataLoader(test_dataset, batch_size=1, follow_batch=batch_vars)


# Baseline optimizer:
optimizer = torch.optim.Adam(net.parameters(), lr=3e-4, amsgrad=True)
best_loss = 1e10  # We save the "best model so far"

starting_epoch = 0
if args.restart_training != "":
    #checkpoint = torch.load("models/" + args.restart_training)
    checkpoint = torch.load(args.restart_training)
    net.load_state_dict(checkpoint["model_state_dict"])
    optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
    starting_epoch = checkpoint["epoch"]
    best_loss = checkpoint["best_loss"]

early_stopping = EarlyStopping(patience=args.early_stopping, verbose=True)

# Training loop (~100 times) over the dataset:
for i in range(starting_epoch, args.n_epochs):
    # Train first, Test second:
    for dataset_type in ["Train", "Validation", "Test"]:
        if dataset_type == "Train":
            test = False
        else:
            test = True

        suffix = dataset_type
        if dataset_type == "Train":
            dataloader = train_loader
        elif dataset_type == "Validation":
            dataloader = val_loader
        elif dataset_type == "Test":
            dataloader = test_loader

        # Perform one pass through the data:
        info = iterate(
            net,
            dataloader,
            optimizer,
            args,
            test=test,
            summary_writer=writer,
            epoch_number=i,
            flex = args.flexibility,
        )

        if dataset_type == "Validation" and i > 49:
            early_stopping(info["Loss"])

            if early_stopping.early_stop:
                print("Early stopping")
                exit()

        # Write down the results using a TensorBoard writer:
        for key, val in info.items():
            if key in [
                "Loss",
                "ROC-AUC",
                "Distance/Positives",
                "Distance/Negatives",
                "Matching ROC-AUC",
            ]:
                writer.add_scalar(f"{key}/{suffix}", np.mean(val), i)

            if "R_values/" in key:
                val = np.array(val)
                writer.add_scalar(f"{key}/{suffix}", np.mean(val[val > 0]), i)

        if dataset_type == "Validation":  # Store validation loss for saving the model
            val_loss = np.mean(info["Loss"])

    if True:  # Additional saves
        if val_loss < best_loss:
            print("Validation loss {}, saving model".format(val_loss))
            torch.save(
                {
                    "epoch": i,
                    "model_state_dict": net.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "best_loss": best_loss,
                },
                model_path + "_epoch{}".format(i),
            )

            best_loss = val_loss