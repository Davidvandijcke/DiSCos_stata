program define gen_data
    version 16.0
	
	clear
	
    * If you want the exact same data each time:
    set seed 12345

    //------------------------------------------------------
    // Create artificial dataset with varying distributions
    //------------------------------------------------------
    * Step 1: Generate IDs and Time Periods
    set obs 20                         // Number of IDs
    gen id = _n                        // Create unique IDs
    expand 20                          // Duplicate each ID 20 times (for time periods)
    bysort id: gen time = _n          // Generate time within each ID
    expand 50                          // Create 50 observations per ID-time pair

    * Step 2: Generate group-specific means and standard deviations
    bysort id time: gen double group_mean = runiform()*10 - 5   
    bysort id time: gen double group_sd = runiform()*2 + 0.5    

    * Step 3: Generate the y variable with group-specific means and variances
    gen double y = group_mean + group_sd * rnormal()

    // Add treatment effect that varies across the distribution
    replace y = y + 1 + 0.5*y if id==1 & time>=10  // Treatment effect increases with y
end
