# Effects of Smoking Analysis

## Overview
This repository contains the code and data used for analyzing gene expression data to identify genes that respond differently to smoking in men vs. women.

## Features
- Conducts 2-way ANOVA to analyze gene expression data.
- Visualizes p-values using a histogram.
- Applies False Discovery Rate (FDR) correction for multiple testing.
- Identifies significant genes and intersects them with gene lists related to specific biological processes.
## Getting Started

### Prerequisites

- [Python](https://www.python.org/) and [pip](https://pip.pypa.io/)
- [TensorFlow](https://www.tensorflow.org/) and [Keras](https://keras.io/)
- [NumPy](https://numpy.org/)

### Installation

1. Clone the repository:

   git clone [Effects-of-Smoking](https://github.com/ugendar07/Effects-of-Smoking.git)

## Usage
- Prepare your dataset for model, ensuring that it includes samples.

- Configure the model parameters and dataset paths in the configuration files.

- Train the model
  
- Evaluate the trained model


## Data Set
The input data consists of gene expression data from both male and female individuals, categorized by smoking status. Additionally, gene lists related to Xenobiotic Metabolism, Free Radical Response, DNA Repair, and Natural Killer Cell Cytotoxicity are provided for comparison.
- [Xenobiotic Metabolism](https://drive.google.com/open?id=1hRAKYSvN6mNa4DWZxmLNP5PnjqFuGZtJ)
- [Free Radical Response](https://drive.google.com/open?id=16Ot-Kgmyvs-yNBDhobRsy5DSWNvZXu6y)
- [DNA Repair](https://drive.google.com/file/d/16b5kBgvLmCSilzfNTD-p2EOTZhpxbnSV/view?usp=sharing)
- [Natural Killer Cell Cytotoxicity](https://drive.google.com/file/d/1vykPnfqafHkSd1ivxBXaORjvPkfZ3IE4/view?usp=sharing)


## Acknowledgments
I acknowledge the use of the 2-way ANOVA framework for analyzing gene expression data in this project.


## Contact
For questions or inquiries, please contact [ugendar](mailto:ugendar07@gmail.com) .
