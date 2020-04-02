class SynFormula:
    """
    data strcuture to store information for a synthetic species
    """

    def __init__(self, formula, covariate_list, coefficient_dict, new_values, probability):
        self.formula = formula
        self.covariates = covariate_list
        self.coefficients = coefficient_dict
        self.gen_values = new_values
        self.probability = probability
