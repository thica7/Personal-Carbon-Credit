# Step 1
def calculate_carbon_cap(cur_carbon,T,P):
    carbon_cap = cur_carbon/T/P
    return carbon_cap


# Step 2
def calculate_carbon_baseline(carbon_cap):
    carbon_baseline = 0.9*carbon_cap
    return carbon_baselin

# Step 3
def calculate_carbon_allocation(carbon_baseline):
    carbon_allocation = beta * carbon_baseline
    # Where beta is based on specific situations.
    return carbon_allocation
