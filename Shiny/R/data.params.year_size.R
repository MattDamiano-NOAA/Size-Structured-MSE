# Dolphinfish MSE
# Dolphinfish OM: Input data file
# Last modified: 10-11-2023 by M. Damiano

initial_year <- 1986
terminal_year <- 2022

years = rep(initial_year:terminal_year) # set the number of years
years_count = length(years) # set annual time step
year_start = years[1] # set start year
year_end = tail(years, n = 1) # a bit redundant, but there it is

seasons_count = 4 # set seasonal timestep; 1=no seasonality
n_block = 1 # set blocks for time-varying quantities, parameters; 1=no blocks

regions_count = 7 # number of regions for spatially-explicit model (we may want 8 to split the Caribbean into north and south later)
size_bins = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000) # define size bins
size_bm = c(150, 250, 350, 450, 550, 650, 750, 850, 950, 1050, 1150, 1250, 1350, 1450, 1550, 1650, 1750, 1850, 1950, 2050) # midpoints
size_bins_count = length(size_bm) # set length structure variable based on number of size bins

# What fisheries do we know we need coastwide?
# Commercial (maybe combined or just longline is ok);
# Private recreational x 2
# For-hire recreational x 2
# High seas intl (eventually)
# Discard fleet
# Then we will probably need a Venezuelan artisanal drift gill net for lower Caribbean
F_init_v <- c(0.1, 0.3, 0.3, 0.05, 0.05, 0.05)
n_fleets = length(F_init_v) # Fleets: commercial, private rec N, private rec S, for-hire N, for-hire S, discard
