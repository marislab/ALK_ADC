# ALK Expression Outlier Analysis

This code identifies histologies with outlier ALK expression and generates figures showing ALK expression in those histologies, along with expression in normal tissues from GTEx. The plot generated by this script corresponds to Figure 1 A & B.

## Usage

To run the script:

```sh
./generateFigure1AB.R
```

The output will be a plot saved in your project working directory with the other plots generated by the script.

**Note:**

* Make sure you have read/write permissions to the folder where you execute the script.
* The data should be present in the same folder as the script, as it is designed to read data from the current directory.
* Output plots will be saved in your current working directory.

