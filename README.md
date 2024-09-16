# Improving Antibody-Antigen Interaction Prediction through Flexibility with ESMFold #

Antibodies are essential proteins in the immune system due to their capacity to bind to specific antigens. They also play a critical role in developing vaccines and treatments for infectious diseases. Their complex structure, with variable regions for antigen binding and flexible hinge regions, presents challenges for accurate computational modeling. Recent advancements in deep learning have revolutionized protein structure prediction. Despite these advancements, predicting interactions between antibodies and antigens remains challenging, mainly due to the flexibility of antibodies and the dynamic nature of binding events. In this study, we employ fingerprint-based methodologies, incorporating ESMFold confidence scores as a flexibility feature, to model the flexibility of Ab-Ag interactions.With our methodology, we show how including flexibility has improved Ab-Ag interactions by 0.8% arriving at an AUC-ROC of 90%.


## Environment ##

To create the environment, run the following command:

'''
conda create -n flexibility python=3.8.18 pip -y
conda activate flexibility
pip install -r requirements.txt
'''

## Model training and evaluation ##
To train the model and to evaluate the model run the following commands:

'''
python main_training.py
python main_inference.py
'''
