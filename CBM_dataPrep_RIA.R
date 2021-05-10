defineModule(sim, list(
  name = "CBM_dataPrep_RIA",
  description = "A data preparation module to format and prepare user-provided input to the SpaDES forest-carbon modelling familly.",
  keywords = NA,
  authors = c(
    person("Celine", "Boisvenue", email = "Celine.Boisvenue@canada.ca", role = c("aut", "cre"))
  ),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.2", CBM_dataPrep_RIA = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "CBM_dataPrep_RIA.Rmd"),
  reqdPkgs = list(
    "data.table", "fasterize", "magrittr", "raster", "RSQLite", "sf",
    "PredictiveEcology/CBMutils (>= 0.0.6)",
    "PredictiveEcology/LandR@development"
  ),
  parameters = rbind(
    # defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(
      ".plotInitialTime", "numeric", NA, NA, NA,
      "This describes the simulation time at which the first plot event should occur"
    ),
    defineParameter(
      ".plotInterval", "numeric", NA, NA, NA,
      "This describes the simulation time interval between plot events"
    ),
    defineParameter(
      ".saveInitialTime", "numeric", NA, NA, NA,
      "This describes the simulation time at which the first save event should occur"
    ),
    defineParameter(
      ".saveInterval", "numeric", NA, NA, NA,
      "This describes the simulation time interval between save events"
    ),
    defineParameter(
      ".useCache", "logical", FALSE, NA, NA,
      paste(
        "Should this entire module be run with caching activated?",
        "This is generally intended for data-type modules,",
        "where stochasticity and time are not relevant"
      )
    )
  ),
  inputObjects = bindrows(
    expectsInput(
      objectName = "cbmData", objectClass = "dataset",
      desc = "S4 object created from selective reading in of cbm_default.db in CBM_defaults module",
      sourceURL = NA
    ),
    expectsInput(
      objectName = "pooldef", objectClass = "character",
      desc = "Vector of names (characters) for each of the carbon pools, with `Input` being the first one",
      sourceURL = NA
    ),
    expectsInput(
      objectName = "PoolCount", objectClass = "numeric",
      desc = "count of the length of the Vector of names (characters) for each of the carbon pools, with `Input` being the first one",
      sourceURL = NA
    ),
    expectsInput(objectName = "dbPath", objectClass = "character", desc = NA, sourceURL = NA),
    expectsInput(objectName = "sqlDir", objectClass = "character", desc = NA, sourceURL = NA),
   expectsInput(
      objectName = "cbmAdmin", objectClass = "dataframe",
      desc = "Provides equivalent between provincial boundaries, CBM-id for provincial boundaries and CBM-spatial unit ids",
      sourceURL = "https://drive.google.com/file/d/1xdQt9JB5KRIw72uaN5m3iOk8e34t9dyz"
    ),
    expectsInput(
      objectName = "userDistFile", objectClass = "character",
      desc = paste("User provided file name that identifies disturbances for simulation",
                   "(key words for searching CBM files, if not there the userDist will be created with defaults"),
      sourceURL = NA
    ),
    expectsInput(
      objectName = "userDist", objectClass = "data.table",
      desc = "User provided file that identifies disturbances for simulation (distName),
      raster Id if applicable, and wholeStand toggle (1 = whole stand disturbance, 0 = partial disturbance),
      if not there it will use userDistFile",
      sourceURL = "https://drive.google.com/file/d/1Gr_oIfxR11G1ahynZ5LhjVekOIr2uH8X"
    ),
        expectsInput(
      objectName = "userGcM3File", objectClass = "character",
      desc = paste("User-provided pointer to the file containing: GrowthCurveComponentID,Age,MerchVolume.",
                   "Default name userGcM3"),
      sourceURL = NA
    ),
    expectsInput(
      objectName = "userGcM3", objectClass = "dataframe",
      desc = "User file containing: GrowthCurveComponentID,Age,MerchVolume. Default name userGcM3",
      sourceURL = "https://drive.google.com/file/d/1BYHhuuhSGIILV1gmoo9sNjAfMaxs7qAj"
    ),
    expectsInput(
      objectName = "masterRaster", objectClass = "raster",
      desc = "Raster built in based on user provided info. Will be used as the raster to match for all operations"
    ),
    expectsInput(
      objectName = "allPixDT", objectClass = "data.table",
      desc = "Data table built for all pixels (incluing NAs) for the four essential raster-based information,
      growth curve location (gcID), ages, ecozones and spatial unit id (CBM-parameter link)"
    ),
    expectsInput(
      ## URL RIA FOR scfmFires rasters is this https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8
      ## URL below is for the data table for the 526 years of scfm fire sims
      objectName = "disturbanceRasters", objectClass = "dataframe",
      desc = "RIA 2020 specific - fires rasters were too big forlow RAM machines. Created a data table for with pixel burnt and year of burn",
      sourceURL = "https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M"
    ),
   expectsInput(
     objectName = "distIndexDT", objectClass = "data.table",
     desc = "Data table built in case the disturbanceRaster data.table was built on a different raster then the one we use for simulations"
    )
   ),
  outputObjects = bindrows(
    createsOutput(objectName = "pools", objectClass = "matrix", desc = NA),
    createsOutput(objectName = "curveID", objectClass = "character",
                  desc = "Vector of column names that together, uniquely define growth curve id"),
    createsOutput(
      objectName = "ages", objectClass = "numeric",
      desc = "Ages of the stands from the inventory in 1990"
    ),
    createsOutput(
      objectName = "nStands", objectClass = "numeric",
      desc = "not really the number of stands, but the number of pixel groups"
    ),
    createsOutput(
      objectName = "gcids", objectClass = "numeric",
      desc = "The identification of which growth curves to use on the specific stands provided by..."
    ),
    createsOutput(
      objectName = "historicDMIDs", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating historical disturbance type, linked to the S4 table called cbmData. Only Spinup."
    ),
    createsOutput(
      objectName = "lastPassDMIDS", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating final disturbance type, linked to the S4 table called cbmData. Only Spinup."
    ),
    createsOutput(
      objectName = "delays", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating regeneration delay post disturbance. Only Spinup."
    ),
    createsOutput(
      objectName = "minRotations", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating minimum number of rotations. Only Spinup."
    ),
    createsOutput(
      objectName = "maxRotations", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating maximum number of rotations. Only Spinup."
    ),
    createsOutput(
      objectName = "returnIntervals", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating the fixed fire return interval. Only Spinup."
    ),
    createsOutput(
      objectName = "spatialUnits", objectClass = "numeric",
      desc = "The id given to the intersection of province and ecozones across Canada, linked to the S4 table called cbmData"
    ),
    createsOutput(
      objectName = "ecozones", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating the numeric represenation of the Canadian ecozones, as used in CBM-CFS3"
    ),
    createsOutput(
      objectName = "level3DT", objectClass = "data.table",
      desc = paste("the table linking the spu id, with the disturbance_matrix_id and the events.",
                   "The events are the possible raster values from the disturbance rasters of Wulder and White.")
    ),
    createsOutput(
      objectName = "spatialDT", objectClass = "data.table",
      desc = "the table containing one line per pixel"
    )
  )
))

doEvent.CBM_dataPrep_RIA <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "CBM_dataPrep_RIA", "save")
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "CBM_dataPrep_RIA", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
      "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
      sep = ""
    ))
  )
  return(invisible(sim))
}


Init <- function(sim) {
  ## making sure the CBM_defaults.R module was run
  io <- inputObjects(sim, currentModule(sim))
  objectNamesExpected <- io$objectName
  available <- objectNamesExpected %in% ls(sim)

  omit <- which(objectNamesExpected %in% c("userDistFile", "userGcM3File"))
  available <- available[-omit]

  if (any(!available)) {
    stop(
      "The inputObjects for CBM_dataPrep are not all available:",
      "These are missing:", paste(objectNamesExpected[!available], collapse = ", "),
      ". \n\nHave you run ",
      paste0("CBM_", c("defaults"), collapse = ", "),
      "?"
    )
  }

  ### RIA 2020: we created a dt instead of a series of rasters in the .inputObjects.
#
#   age <- sim$ageRaster
#   gcIndex <- sim$gcIndexRaster
#   spuRaster <- sim$spuRaster # made in the .inputObjects
#   ecoRaster <- sim$ecoRaster # made in the .inputObjects
#   ## End rasters------------------------------------------------------------------


  ## Create the data table of all pixels and all values for the study area----------------
  # level2DT <- data.table(
  #   spatial_unit_id = spuRaster[], ages = age[], pixelIndex = 1:ncell(age),
  #   growth_curve_component_id = gcIndex[], growth_curve_id = gcIndex[],
  #   ecozones = ecoRaster[]
  # )
  # keep only the pixels that have all the information: the pixels that will be simulated
  spatialDT <- sim$allPixDT[!is.na(ages),]

  #spatialDT <- level2DT
  ## END data.table of all pixels---------------------------------------------------------


  ## Create the pixel groups: groups of pixels with the same attributes ---------------
  setkeyv(spatialDT, "pixelIndex")
  spatialDT$pixelGroup <- Cache(LandR::generatePixelGroups, spatialDT,
    maxPixelGroup = 0,
    columns = c("ages", "spatial_unit_id", "growth_curve_component_id", "ecozones")
  )
  setkeyv(spatialDT, "pixelIndex")

  spatialDT <- spatialDT[, .(
    ages, spatial_unit_id, pixelIndex,
    growth_curve_component_id, growth_curve_id, ecozones, pixelGroup
  )]
  setkeyv(spatialDT, "pixelIndex")
  # spatialDT <- spatialDT[order(pixelIndex), ]
  sim$spatialDT <- spatialDT
  # end create pixel groups-------------

  ## Data.table for simulations (one row per pixel group)---------------------
  # this table will be the pixel groups that are used in the spinup procedure in
  # the CBM_core spinup event

  level3DT <- unique(spatialDT[, -("pixelIndex")])
  setkeyv(level3DT, "pixelGroup")

  sim$curveID <- c("growth_curve_component_id", "ecozones") # "id_ecozone" # TODO: add to metadata -- use in multiple modules
  curveID <- sim$curveID
  sim$gcids <- factor(gcidsCreate(level3DT[, ..curveID]))
  set(level3DT, NULL, "gcids", sim$gcids)
  sim$level3DT <- level3DT
  ## End data.table for simulations-------------------------------------------


  ## TODO: problem with ages<=1
  ##################################################### # SK example: can't seem
  #to solve why growth curve id 52 (white birch, good # productivity) will not
  #run with ages= c(0,1,2) it gets stuck in the spinup. Tried ages==1, # and
  #ages==2. Maybe because the first few years of growth are 0 ? (to check) it #
  #does not grow and it does not fill-up the soil pools. # Notes: the GAMs are
  #fit on the cumulative curves of carbon/ha for three # pools. This is to make
  #sure the curves go through 0...but maybe it would # work better for GAMs to
  #be fit on the increments (?). # since all growth curves are for merchantible
  #timber (with diameter limits), it is acceptable to start all increments at
  #the level of year==3.
  #work for this problem for most curves for now: this is from SK runs
  #sim$level3DT[ages==0 & growth_curve_component_id==52,ages:=3]
 ######################################
  ##################### temp fix should

  #sim$level3DT[ages <= 1, ages := 3]

  setorderv(sim$level3DT, "pixelGroup")

  ## Creating all the vectors for the spinup --------------------------------
  sim$ages <- sim$level3DT[, ages]
  sim$nStands <- length(sim$ages)
  sim$pools <- matrix(ncol = sim$PoolCount, nrow = sim$nStands, data = 0)
  colnames(sim$pools) <- sim$pooldef
  sim$pools[, "Input"] <- rep(1.0, nrow(sim$pools))
  #curveID <- sim$curveID
  #sim$gcids <- as.integer(sim$level3DT[, ..curveID][[sim$curveID]])
  sim$delays <- rep.int(0, sim$nStands)
  sim$minRotations <- rep.int(10, sim$nStands)
  sim$maxRotations <- rep.int(1000, sim$nStands)
  setkeyv(sim$level3DT, "spatial_unit_id")
  spinupParameters <- as.data.table(sim$cbmData@spinupParameters[, c(1, 2)])
  setkeyv(spinupParameters,"spatial_unit_id")
  retInt <- merge.data.table(sim$level3DT, spinupParameters,
                  by = "spatial_unit_id", all.x = TRUE)

  setkeyv(retInt, "pixelGroup")
  setkeyv(sim$level3DT, "pixelGroup")
  sim$returnIntervals <- retInt[, "return_interval"]
  sim$spatialUnits <- sim$level3DT[, spatial_unit_id]
  sim$ecozones <- sim$level3DT$ecozones

  ################################################################################
  ## matching the disturbances with the Disturbance Matrix IDs in CBM-CFS3 defaults
  ################################################################################
  # Matching disturbances to CBM disturbance matrix id---------------------------------

  ## WBI RIA: skipping the SK building of mySpuDmids.
  ## In this case the user provided complete userDist.


                # # make the disturbance look-up table to the disturbance_matrix_id(s)
                # # making sim$mySpuDmids
                #
                # userDist <- sim$userDist
                #
                # # Most cases will only require fire (wildfire) and a clearcut. There are 426
                # # disturbance matrices identified in the archive of CBM
                # # (sim$cbmData@disturbanceMatrix). Matrices are associated with spatial units
                # # (sim$cbmData@disturbanceMatrixAssociation). User can select any disturbance
                # # they want to represent. Some disturbance matrices are based on data but most
                # # are expert opinion in the CBM-CFS3 archive.
                # # Disturbance Matrices are specific to spatial_units_ids--------------
                # spu <- unique(sim$spatialDT$spatial_unit_id)
                # # what disturbances in those spu(s)?
                # # spuDist() function is in CBMutils package
                # # it lists all the possible disturbances in the CBM-CFS3 archive for that/those
                # # spatial unit with the name of the disturbance in the 3rd colum.
                # listDist <- spuDist(spu, sim$dbPath)
                #
                #
                # ## Example specific for SK (as per Boisvenue et al 2016)
                # # Disturbances are from White and Wulder and provided as yearly rasters
                # # raster values 1 to 5
                # # #C:\Celine\GitHub\CBM_\data\forIan\SK_data\disturbance_Sask\ReadMe.txt
                # # # Fire =  1
                # # # Harvest = 2
                # # # Lcondition = 3
                # # # Road = 4
                # # # Unclass = 5
                # # Whatever number of disturbances identified that will be used in the
                # # simulation, each disturbance has to have one one disturbance matrix id
                # # associated with it.
                # # make mySpuDmids (distNames,rasterId,spatial_unit_id,disturbance_matrix_id)
                # distName <- c(rep(userDist$distName, length(spu)))
                # rasterId <- c(rep(userDist$rasterId, length(spu)))
                # wholeStand <- c(rep(userDist$wholeStand, length(spu)))
                # spatial_unit_id <- c(sort(rep(spu, length(userDist$distName))))
                # mySpuDmids <- data.table(distName, rasterId, spatial_unit_id, wholeStand)
                #
                # dmid <- data.frame(spatial_unit_id = integer(), disturbance_matrix_id = integer())
                #
                # for (i in 1:length(mySpuDmids$distName)) {
                #   ### DANGER HARD CODED FIXES
                #   ## to do: present the user with options that live in listDist for the
                #   ## specific spu or in sim$cbmData@disturbanceMatrix
                #   if (mySpuDmids$distName[i] == "clearcut") {
                #     dmid[i, ] <- cbind(mySpuDmids$spatial_unit_id[i], 409)
                #   } else {
                #     getDist <- listDist[grep(mySpuDmids$distName[i], listDist[, 3], ignore.case = TRUE), 1:2]
                #     getDist <- getDist[getDist$spatial_unit_id == mySpuDmids$spatial_unit_id[i], ]
                #     dmid[i, ] <- getDist[1, ]
                #   }
                # }
                #
                # ## bunch of warnings here...
                # mySpuDmids <- data.table(mySpuDmids, dmid$disturbance_matrix_id)
                # names(mySpuDmids) <- c("distName", "rasterId", "spatial_unit_id", "wholeStand", "disturbance_matrix_id")
                # sim$mySpuDmids <- mySpuDmids
  # sim$mySpuDmids <- sim$userDist

  # need to match the historic and last past dist to the spatial unit
  # DECISION: both the last pass and the historic disturbance will be the same
  # for the fire return interval runs

  ## TO DO: in Canada historic DMIDs will always be fire, but the last past may
  ## not, it could be harvest. Make this optional and give the user a message
  ## saying these are the defaults.


  mySpuFires <- sim$userDist[grep("wildfire", sim$userDist$distName, ignore.case = TRUE), ]

  myFires <- mySpuFires[spatial_unit_id %in% unique(sim$level3DT$spatial_unit_id), ]
  setkey(myFires, spatial_unit_id)
  setkey(sim$level3DT, spatial_unit_id)
  # this is mainly to make them the same length at the number of pixel groups
  histLastDMIDs <- merge(sim$level3DT, myFires)
  sim$historicDMIDs <- histLastDMIDs$disturbance_matrix_id
  ## TO DO: this is where it could be something else then fire
  sim$lastPassDMIDS <- histLastDMIDs$disturbance_matrix_id

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # CBM defaults ----------------------------------------------------------
  # if we chose to not use the RSQLite library in this module, and extract
  # disturbance matrix id (dmid) from sim$cbmData@disturbanceMatrixAssociation,
  # then $sqlDir and $dbPath are not needed.
  if (!suppliedElsewhere(sim$sqlDir)) {
    sim$sqlDir <- file.path(dPath, "cbm_defaults")
  }
  if (!suppliedElsewhere(sim$dbPath)) {
    sim$dbPath <- file.path(dPath, "cbm_defaults", "cbm_defaults.db")
  }

  if (!suppliedElsewhere(sim$cbmData)) {
    spatialUnitIds <- as.matrix(getTable("spatialUnitIds.sql", sim$dbPath, sim$sqlDir))
    disturbanceMatrix <- as.matrix(getTable("disturbanceMatrix.sql", sim$dbPath, sim$sqlDir))
    sim$cbmData <- new("dataset",
      turnoverRates = as.matrix(getTable("turnoverRates.sql", sim$dbPath, sim$sqlDir)),
      rootParameters = as.matrix(getTable("rootParameters.sql", sim$dbPath, sim$sqlDir)),
      decayParameters = as.matrix(getTable("decayParameters.sql", sim$dbPath, sim$sqlDir)),
      spinupParameters = as.matrix(getTable("spinupParameters.sql", sim$dbPath, sim$sqlDir)),
      climate = as.matrix(getTable("climate.sql", sim$dbPath, sim$sqlDir)),
      spatialUnitIds = spatialUnitIds,
      slowAGtoBGTransferRate = as.matrix(0.006),
      biomassToCarbonRate = as.matrix(0.5),
      stumpParameters = as.matrix(getTable("stumpParameters.sql", sim$dbPath, sim$sqlDir)),
      overmatureDeclineParameters = as.matrix(getTable("overmaturedecline.sql", sim$dbPath, sim$sqlDir)),
      disturbanceMatrix = disturbanceMatrix,
      disturbanceMatrixAssociation = as.matrix(getTable("disturbanceMatrixAssociation.sql", sim$dbPath, sim$sqlDir)),
      disturbanceMatrixValues = as.matrix(getTable("disturbanceMatrixValues.sql", sim$dbPath, sim$sqlDir)),
      landclasses = as.matrix(getTable("landclasses.sql", sim$dbPath, sim$sqlDir)),
      pools = as.matrix(getTable("pools.sql", sim$dbPath, sim$sqlDir)),
      domPools = as.matrix(getTable("domPools.sql", sim$dbPath, sim$sqlDir))
    )
  }
  if (!suppliedElsewhere(sim$pooldef)) {
    sim$pooldef <- CBMutils::.pooldef
    sim$PoolCount <- length(sim$pooldef)
  }
  ## END CBM defaults-----------------------------------------------------

  # user provided data tables (3)------------------------------------------------------

  # 1. growth and yield information
  # userGcM3 and userGcM3File, these files are the m3/ha and age info by growth
  # curve ID, columns should be GrowthCurveComponentID	Age	MerchVolume
  ## TODO add a data manipulation to adjust if the m3 are not given on a yearly basis
  if (!suppliedElsewhere("userGcM3", sim)) {

    if (!suppliedElsewhere("userGcM3File", sim)) {
      sim$userGcM3File <- extractURL("userGcM3")
    }

    sim$userGcM3 <- prepInputs(url = sim$userGcM3File,
                               fun = "data.table::fread",
                               destinationPath = dPath,
                               #purge = 7,
                               filename2 = "curve_points_table.csv")

    # message(
    #   "User has not supplied growth curves (m3 by age or the file name for the growth curves). ",
    #   "The default will be used which is for a region in Saskatchewan."
    # )
    ## RIA 2020 specific
    sim$userGcM3[, V1 := NULL]
    names(sim$userGcM3) <- c("GrowthCurveComponentID", "Age", "MerchVolume")
  }

  # 2. Disturbance information - see disturbance raster below
  # this may be provided by the user, by the defaults or by other modules/family
  # of modules. It is the link between the spatial location of the disturbance
  # (like a raster value) and the disturbance name.
  if (!suppliedElsewhere(sim$userDist)) {
    if (!suppliedElsewhere(sim$userDistFile)) {
      sim$userDistFile <- extractURL("userDist")
    }
    sim$userDist <- prepInputs(url = sim$userDistFile,
                                 fun = "data.table::fread",
                                 destinationPath = dPath,
                                 # purge = 7,
                                 filename2 = "mySpuDmids.csv")
  }

  # 3. cbmAdmin needed to create spuRaster below (a bit convoluted ## TODO )
  if (!suppliedElsewhere("cbmAdmin", sim)) {
    sim$cbmAdmin <- prepInputs(url = extractURL("cbmAdmin"),
                               fun = "data.table::fread",
                               destinationPath = dPath,
                               #purge = 7,
                               filename2 = "cbmAdmin.csv")
  }

  # END user provided data tables (3)------------------------------------------------------

  # user provided rasters or spatial information----------------------------------------

  ## Rasters
  ## user provides raster to match (masterRaster) which is a raster for the
  ## study area, it will define the crs etc, for all other layers. The user also
  ## provides age raster, and a raster linking each growth curve to pixels (gcIndex).
  ## Using the masterRaster, the ecozone raster is made (Canadian ecozones) and the
  ## spatial unit raster. The spatial units are a CBM-CFS3 specific location
  ## that is the intersection of the ecozones and administrative boundaries.
  ## These spatial units (or spu) and the ecozones link the CBM-CFS3 ecological
  ## parameters to the right location (example: decomposition rates).
  ##

  #1. VRI rasters gcID and age
  # study area raster
  RIArtm <- prepInputs(url = "https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing",
                       destinationPath = dPath)
  # This works for the fire return interval runs
  #forest inventory info. this makes a raster stack of the two main rasters, gcIDRaster and ageRaster
  RIA_VRIstack <- Cache(prepInputsVRI,VRIurl = "https://drive.google.com/file/d/1LXSX8M46EnsTCM3wGhkiMgqWcqTubC12",
                        dPath = dPath,
                        rasterToMatch = RIArtm
  )
  names(RIA_VRIstack) <- c("gcIDRaster", "ageRaster")


  #2. Seperate the rasters in the stack
  gcIDRaster <- RIA_VRIstack$gcIDRaster
  ageRaster <- RIA_VRIstack$ageRaster # this is VRI2020 age raster


  #3. Make a masterRaster and make sure there are no NAs
  masterRaster <- raster::raster(gcIDRaster)
  ageWData <- !is.na(ageRaster[])
  gcIDWData <- !is.na(gcIDRaster[])
  masterRaster[ageWData] <- gcIDRaster[ageWData]
  gcIDRaster[!ageWData] <- NA


  ## need to check if the VRI2015 age has values in all the same pixels. It
  ## would be better to use the same masterRaster for all three RIA sims sets,
  ## fireReturnInterval,presentDay and the two harvest scenarios.
  ## TODO: the prespInputsVRI function makes a stack of two rasters but in the
  ## case of the presentDay sims, we need a different age (DONE, created the
  ## function prepInputVIRage). This age raster starts from the ESRI file here
  ## https://pub.data.gov.bc.ca/datasets/02dba161-fdb7-48ae-a4bb-bd6ef017c36d/2015/VEG_COMP_LYR_L1_POLY_2015.gdb.zip
  ## which is the 2015VRI for the province of BC.

### changing th order so we can build the PresentDat ageRaster1985
  #4. Make the ecozone Raster (ecoRaster)"http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"
  ecozone <- Cache(prepInputsEcozones, url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
                   dPath = dPath,
                   masterRaster = masterRaster)

  #5. Get just BC
  #provs <- getData("GADM", country = "CAN", level = 1)

  #bc <- provs[provs$NAME_1 == "British Columbia",]

  #6. mathc the SPU values to the ecozone values

  cbmAdminThisSA <- sim$cbmAdmin[adminName == "British Columbia", ]

  rows <- match(ecozone[], cbmAdminThisSA$EcoBoundaryID)
  spatialUnitID <- cbmAdminThisSA[rows,"SpatialUnitID"]

  dtRasters <- as.data.table(cbind(growth_curve_component_id = gcIDRaster[], ages = ageRaster[], ecozones = ecozone[], spatialUnitID))

  # assertion -- if there are both NAs or both have data, then the colums with be the same, so sum is either 0 or 2
  if (isTRUE(P(sim)$doAssertions)) {
    bbb <- apply(dtRasters, 1, function(x) sum(is.na(x)))
    if (!all(names(table(bbb)) %in% c("0", "4")))
      stop("should be only 0 or 4s")
  }

  #7. Passed assertion? Save the masterRaster and data.table with all pixel indices.
  sim$masterRaster <- masterRaster

  sim$allPixDT <- as.data.table(cbind(dtRasters, pixelIndex = 1:ncell(gcIDRaster), growth_curve_id = gcIDRaster[]))
  setnames(sim$allPixDT,"SpatialUnitID", "spatial_unit_id")

  # 8. Disturbance rasters. The default example (SK) is a list of rasters, one for
  # # each year. But these can be provided by another family of modules in the
  # # annual event.
  ### TODO give options to the user to provide a raster a data table, a raster list or a raster stack

  if (!suppliedElsewhere("disturbanceRasters", sim)) {
    ## this case is reading in a sparseDT.
    # RTM that the datatable was created with
    RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                          destinationPath = dPath) #you probably already have this raster - RIA_RTM.tif on google
    # need the masterRaster (sim$masterRaster)
    # sparseDT
    # scfmAnnualBurns <- prepInputs(url = 'https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M/view?usp=sharing',
    #                               destinationPath = 'inputs',
    #                               overwrite = TRUE,
    #                               fun = 'readRDS')
    presentDayBurns <- prepInputs(url = 'https://drive.google.com/file/d/1MjQ5y9Txr1ezv0BatG4b_6FpM6wti1b5/view?usp=sharing',
                                  destinationPath = 'inputs',
                                  overwrite = TRUE,
                                  fun = 'readRDS')
    presentDayHarvest <- prepInputs(url = 'https://drive.google.com/file/d/1Ca-kPun7_VF2xV8s40IJ9hri3Ok-6qx5/view?usp=sharing',
                                    destinationPath = 'inputs',
                                    overwrite = TRUE,
                                    fun = 'readRDS')

    # this adjustment is only needed when we are using a different
    IndexRTM <- setValues(RIA_RTM, 1:ncell(RIA_RTM))

    #postProcess to match tempTHLB
    IndexTHLB <- setValues(sim$masterRaster, 1:ncell(sim$masterRaster))

    #postProcess the RTM
    IndexRTM <- postProcess(IndexRTM, rasterToMatch = sim$masterRaster)
    #build matching data.table

    indexDT <- data.table(rtmIndex = getValues(IndexRTM),
                          thlbIndex = getValues(IndexTHLB))
    sim$distIndexDT <- indexDT[!is.na(rtmIndex)]

    # in this case: they are equal.
    # dim(spadesCBMout$distIndexDT)
    # [1] 3112425       2
    # > length(which(spadesCBMout$distIndexDT$rtmIndex == spadesCBMout$distIndexDT$thlbIndex))
    # [1] 3112425

    # presentDay runs will require an extra step: putting the DTs together with
    # a rasterID of 1 for fire and 2 for harvest.
    presentDayBurns[, events := 1L]
    presentDayHarvest[, events := 2L]
    allDist <- rbind(presentDayBurns, presentDayHarvest)
    #the NAs in rtmIndex are pixels that are not in THLB (but inside the landscape) - we can remove them
    sim$disturbanceRasters <- allDist

    ## a function indexAnnualFire() may be used in the annual event
    ## of the CBM_core module to extract the year



    # options(reproducible.useGDAL = FALSE)
    # sim$disturbanceRasters <- Cache(prepInputs,
    #                    url = "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8",
    #                    fun = "raster::brick",
    #                    rasterToMatch = sim$masterRaster,
    #                    datatype = "INT1U",
    #                    useGDAL = FALSE)

    ## TODO this is a brick, with 526 rasters that low RAM system cannot
    ## handle. Instead a sparseDT was create and is read-in here. We need to
    ## add here the capacity to deal with a bunch of rasters in a folder,
    ## (like in the SK runs), a stack of rasters (below), raster brick (above)
    ## or polygons? OR a sparseDT
    # stack
    # sim$disturbanceRasters <- Cache(prepInputs,
    #                               url = "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8",
    #                               fun = "raster::stack",
    #                               rasterToMatch = masterRaster,
    #                               useGDAL = FALSE)
    # rasters in a folder
    # distHere <- extractURL(disturbanceRasters)
    # sim$disturbanceRasters <- list.files(distHere,full.names = TRUE) %>%
    #   grep(., pattern = ".grd$", value = TRUE)
    # # if all fails
    # or just one raster? or polygons?
  }

  # 9 make a new age raster that starts from the VRI2015
  # VRI2015 to back-build the age raster


  # sa <- as(extent(masterRaster), "SpatialPolygons")
  # crs(sa) <- crs(masterRaster)
  # loadAge <- function(x, field = "PROJ_AGE_1") {
  #   # a <- Cache(sf::st_read, x) # I used Cache during my development because this takes 37 minutes to run -- I was sick of running it again and again
  #   a <- sf::st_read(x)
  #   a1 <- a[, field]
  #   return(a1)
  # }
  # a <- Cache(prepInputs,
  #            url = "https://pub.data.gov.bc.ca/datasets/02dba161-fdb7-48ae-a4bb-bd6ef017c36d/2015/VEG_COMP_LYR_L1_POLY_2015.gdb.zip",
  #            #fun = quote(loadAge(x = targetFilePath,
  #            #field = "PROJ_AGE_1")),
  #            targetFile = "VEG_COMP_LYR_L1_POLY_2015.gdb.zip",
  #            archive = NA)
  # #studyArea = sa
  #                 #rasterTomatch = masterRaster)
  # vriAge2015 <- sf::st_read("C:/Celine/github/spadesCBM_RIA/VEG_COMP_LYR_L1_POLY_2015.gdb.zip")
  #
  # b <- sf::st_transform(vriAge2015, st_crs(sa))
  #
  # ageRaster2015 <- fasterize::fasterize(b, masterRaster, field = "PROJ_AGE_1")
  # ## HERE
  # ageRaster2015[] <- as.integer(ageRaster2015[])


  ageRaster2015 <- Cache(prepInputsVRIage,
                   VRIurl = "https://pub.data.gov.bc.ca/datasets/02dba161-fdb7-48ae-a4bb-bd6ef017c36d/2015/VEG_COMP_LYR_L1_POLY_2015.gdb.zip",
                   dPath = dPath,
                   rasterToMatch = masterRaster,
                   targetFile = "VEG_COMP_LYR_L1_POLY_2015.gdb.zip",
                   field = "PROJ_AGE_1")

  age2015 <- ageRaster2015[]
  age2020 <- ageRaster[]
  # NAs 2015 5539822
  # NAs 2020 5435395


  ## figure out what ages the NAs have in 2020.
  ## Note that some ages in the 2020 that have no ages in the 2015 raster AND that
  ## are disturbed (so in the sim$disturbanceRasters data table) are old...older
  ## then the age from the dist...but most are below 50.


  ageDT <- data.table(pixelIndex = sim$allPixDT$pixelIndex, age2015, age2020)
  ageNoMatch <- ageDT[is.na(age2015) & !is.na(age2020),]
  ## NOTEs: there are 104427 more pixels with NAs in age2015 than in age2020.
  ageDT[!is.na(age2015) & is.na(age2020),]
  ## There are no pixels that have ages in 2015 and are NAs-age in 2020.

  # ageNA1985, what is their year of disturbance?

  setkeyv(ageNoMatch, "pixelIndex")
  setkeyv(allDist,"pixelID")

  # do any pixels burn twice?
  length(unique(presentDayBurns$pixelID)) == dim(presentDayBurns)[1]
  #TRUE
  # harvested twice?
  length(unique(presentDayHarvest$pixelID)) == dim(presentDayHarvest)[1]
  #TRUE

  # get only the burnt/harvested pixels in the noMatch
  ageNAburns <- presentDayBurns[pixelID %in% ageNoMatch$pixelIndex, ]
  setnames(ageNAburns,"pixelID", "pixelIndex")
  ageNAharvest <- presentDayHarvest[pixelID %in% ageNoMatch$pixelIndex, ]
  setnames(ageNAharvest,"pixelID", "pixelIndex")
  # Merge them with ageNoMatch
  #dt_a[dt_b, on = .(b = y)]
  ageNA1985 <- merge.data.table(ageNoMatch, ageNAburns, by = "pixelIndex", all.x = TRUE)
  setnames(ageNA1985,"year", "burnYear")
  ageNA1985[, events := NULL]
  ageNA1985 <- merge.data.table(ageNA1985, ageNAharvest, by = "pixelIndex", all.x = TRUE)
  setnames(ageNA1985,"year", "harvestYear")
  ageNA1985[, events := NULL]
  ageNA1985$noDist <- 0
  ageNA1985[which(is.na(burnYear) & is.na(harvestYear)),]$noDist <- 1

  # make a column of "straight substraction"
  ageNA1985[, substract := (age2020 - 35)]

  ## pixels that have no disturbance and are >0 in 1985 get this age ##################
  ageNA1985$PixSubtract1985 <- 999
  ageNA1985[noDist == 1 & substract >= 0,]$PixSubtract1985 <- ageNA1985[noDist == 1 & substract >= 0,]$substract
  ageNA1985[PixSubtract1985 == 999,]
  # still 23177 that are not dealt with.

  ## the ones that have a non-negative age at the time of disturbance will get the
  ## age in the year prior to disturbance
  ageNA1985[, PixFireYrAge := age2020 - (2020-burnYear)]
  ageNA1985[, PixCutYrAge := age2020 - (2020-harvestYear)]
  ## this gives the ages of the burnt pixels in 1985 (values>0)
  ageNA1985[, firePix1985age := PixFireYrAge - (burnYear - 1985)]
  ageNA1985[, cutPix1985age := PixCutYrAge - (harvestYear - 1985)]

  ## how many left?
  ageNA1 <- ageNA1985[PixSubtract1985 != 999] #81250
  ageNA2 <- ageNA1985[ PixSubtract1985 == 999 & (firePix1985age > 0 | cutPix1985age > 0),] #3835

  #104427-85085
  ageNA1985[pixelIndex %in% ageNA1$pixelIndex | pixelIndex %in% ageNA2$pixelIndex] #85085
  probPix <- ageNA1985[!(pixelIndex %in% ageNA1$pixelIndex | pixelIndex %in% ageNA2$pixelIndex)] #19342
  # remove these
  #colsRemove <- c("substract", "PixSubtract1985", "PixFireYrAge")
  #probPix[, c("substract", "PixSubtract1985", "PixFireYrAge") := list(NULL, NULL, NULL)]

  # how many of those are disturbed?
  table(probPix$noDist)
  # 0     1
  # 11183  8159

  # look at disturbed pixels ages in 1985
  # probCuts1985 <- qplot(probPix[noDist == 0,]$cutPix1985age, geom = "histogram")
  # probBurns1985 <- qplot(probPix[noDist == 0,]$firePix1985age, geom = "histogram")

  ## decision: all the disturbed pixels (noDist == 0) that were cut or burnt will
  ## have mature ages at the year of disturbance. Those will be selected from the
  ## age distribution of pixels with the same gcID with ages >80 in age2020.

  distPixToMatch <- probPix[noDist == 0,]
  length(unique(distPixToMatch$pixelIndex)) # no duplicats
  # figure out the gcID for each pixelIndex
  # get all the ages>80 at ages2020 in each gcID
  # randomly select one for each pixelIndex

  #gcIDdistPixNoMatch <- unique(sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex,]$growth_curve_id)
  #84
  rndNumBygcID <- sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex, .(pixelIndex, growth_curve_id)][,.N, by = "growth_curve_id"]

  matchPixDT <- sim$allPixDT[growth_curve_id %in% rndNumBygcID$growth_curve_id & ages>80,.(growth_curve_id, ages)]

  newAgedistPix1 <- matchPixDT[rndNumBygcID, on = "growth_curve_id", nomatch = 0]

  #newPixDT <- newPixDT[, .SD[sample(.N, N)], by = "growth_curve_id"] # take a sample of each growth_curve_id, length N
  #newAgedistPix <- newAgedistPix1[, lapply(.SD, sample(ages, size = N[1])), by = growth_curve_id]

  newAgedistPix <- newAgedistPix1[, .SD[sample(ages,size = N[1])], by = "growth_curve_id"]

  # This is just a test -- does each growth curve id have the same # rows as N says
  newAgedistPix[, .N == N[1],by = "growth_curve_id"]

  # re-attach a pixelIndex
  pixIndDistPixToMatch <- sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex, .(pixelIndex, growth_curve_id)]
  setorder(pixIndDistPixToMatch, growth_curve_id)
  setorder(newAgedistPix, growth_curve_id)
  set(pixIndDistPixToMatch, NULL, "distPixNegAge1985", newAgedistPix$ages)
  pixIndDistPixToMatch[, growth_curve_id := NULL]
  ageNA1985 <- merge.data.table(ageNA1985, pixIndDistPixToMatch, on = "pixelIndex", all.x = TRUE)

  ##
  noDistPixToMatch <- probPix[noDist == 1,]
  ## trying to find the closest pixel with a disturbance in the 1985-2015
  ## disturbance data.table.
  # are there pixels that are disturbed twice in this data?
  countDist <- allDist[, .N, by = pixelID]
  # yes
  # are those in the proPix?
  probPix[pixelIndex %in% countDist[N>1,]$pixelID,] #no

  # create a raster with all the disturbances
  allDistRaster <- raster::raster(masterRaster)
  allDistRaster[] <- NA
  allDistRaster[countDist$pixelID] <- 1

  # trying focal. This takes a raster and gives me back a raster
  f3 <- function(x){
    theNAs <- is.na(x)
    if (all(theNAs))
      NA
    else
      x[sample(which(!theNAs),1)]
  }
  # maybe too small?
  f9 <- focal(allDistRaster,
              w=matrix(1,nrow=3,ncol=3),
              fun = f3)

  # agg that column to the all pixels DT
  set(ageDT, NULL, "f9", f9[])
  # do any of the 3X3 windows cover the last 8159 pixels?
  checkF9 <- ageDT[pixelIndex %in% noDistPixToMatch$pixelIndex,]
  table(checkF9$f9, useNA = "ifany")
  # 1  NaN
  # 4410 3749
  # figure out what year and what disturbances
  # some pixels are disturbed twice but they are not in my probPix
  #
  # make a raster with the dist year
  yrDistRaster <- raster::raster(masterRaster)
  setnames(countDist, "pixelID", "pixelIndex")
  setnames(allDist, "pixelID", "pixelIndex")
  yrDist <- unique(countDist[allDist, on = "pixelIndex", nomatch = 0][,.(pixelIndex, year)])
  yrDistRaster[] <- NA
  yrDistRaster[yrDist$pixelIndex] <- yrDist$year
  yrf9 <- focal(yrDistRaster,
                w=matrix(1,nrow=3,ncol=3),
                fun = f3)
  set(ageDT, NULL, "yrf9", yrf9[])
  # make a raster for the dist type
  eventRaster <- raster::raster(masterRaster)
  eventDist <- unique(countDist[allDist, on = "pixelIndex", nomatch = 0][,.(pixelIndex, events)])
  eventRaster[] <- NA
  eventRaster[eventDist$pixelIndex] <- eventDist$events
  eventf9 <- focal(eventRaster,
                   w=matrix(1,nrow=3,ncol=3),
                   fun = f3)
  set(ageDT, NULL, "eventf9", eventf9[])


  # create a new "distPixels" DT adding the dist to the f9==1
  f9pix <- ageDT[f9==1]
  f9dist <- merge.data.table(allDist, f9pix, all = TRUE)
  # add the growth_curve_id to help the match??
  f9dist <- f9dist[sim$allPixDT, on = 'pixelIndex', nomatch = 0][,.(pixelIndex, year, events,
                                                                    age2015, age2020, f9, yrf9,
                                                                    eventf9, growth_curve_id)]
  f9dist[, targetPix := 0L]
  f9dist[pixelIndex %in% noDistPixToMatch$pixelIndex]$targetPix <- 1

  # are there any pixels with targetPix == 1 and yrf9 !is.na()?
  f9dist[targetPix >0 & !is.na(yrf9)]#4217
  f9FirstHist <- hist(f9dist[targetPix >0 & !is.na(yrf9)]$age2020, plot = FALSE)


  ## create "new" dist to add to the allDist table
  newF9dist <- f9dist[targetPix >0 & !is.na(yrf9), .(pixelIndex, yrf9, eventf9, growth_curve_id)]
  ## assign age at time of dist

  ## calculate age1985

  rndF9BygcID <- newF9dist[,.N, by = "growth_curve_id"]
  #sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex, .(pixelIndex, growth_curve_id)][,.N, by = "growth_curve_id"]

  f9agesDistDT <- sim$allPixDT[growth_curve_id %in% rndF9BygcID$growth_curve_id & ages>80,.(growth_curve_id, ages)]

  newAgeF9atDist <- f9agesDistDT[rndF9BygcID, on = "growth_curve_id", nomatch = 0]

  #newPixDT <- newPixDT[, .SD[sample(.N, N)], by = "growth_curve_id"] # take a sample of each growth_curve_id, length N
  #newAgedistPix <- newAgedistPix1[, lapply(.SD, sample(ages, size = N[1])), by = growth_curve_id]

  newAgeF9atDistPix <- newAgeF9atDist[, .SD[sample(ages,size = N[1])], by = "growth_curve_id"]

  # This is just a test -- does each growth curve id have the same # rows as N says
  newAgeF9atDistPix[, .N == N[1],by = "growth_curve_id"]

  # merge with newF9dist
  setorder(newF9dist, growth_curve_id)
  setorder(newAgeF9atDistPix, growth_curve_id)
  set(newF9dist, NULL, "newAgeF9atDist", newAgeF9atDistPix$ages)
  newF9dist[, ageF91985 := newAgeF9atDist - (yrf9 - 1985)]

  ## add the "new" dist to allDist but make sure there are no NAs
  table(newF9dist$yrf9, useNA = "ifany") # none
  table(newF9dist$eventf9, useNA = "ifany")
  # 1    2  NaN
  # 88 4080   49
  # need to replace those NAs ###HERE
  distVec <- sample(1:2,length(which(is.na(newF9dist$eventf9))), replace = TRUE)
  newF9dist[is.na(eventf9)]$eventf9 <- distVec
  f9AddDist <- newF9dist[,.(pixelIndex, yrf9, eventf9)]
  setnames(f9AddDist,c("yrf9", "eventf9"), c("year", "events"))

  dim(allDist)
  allDist <- rbind(allDist, f9AddDist)

  ## add this age column to the ageNA1985 DT
  f9ages <- newF9dist[,.(pixelIndex, ageF91985)]
  ageNA1985 <- merge.data.table(ageNA1985, f9ages, on = "pixelIndex", all.x = TRUE)

  ## try again with a bigger window
  # what are the pixels left to match?
  toMatch5 <- noDistPixToMatch[!(pixelIndex %in% newF9dist$pixelIndex)]

  ## bigger focal
  # start with yrDistRaster
  yrf25 <- focal(yrDistRaster,
                 w=matrix(1,nrow=7,ncol=7),
                 fun = f3)
  set(ageDT, NULL, "yrf25", yrf25[])
  eventf25 <- focal(eventRaster,
                    w=matrix(1,nrow=7,ncol=7),
                    fun = f3)
  set(ageDT, NULL, "eventf25", eventf25[])

  f25pix <- ageDT[!is.na(yrf25), .(pixelIndex, age2015, age2020, yrf25, eventf25)]
  #f25dist <- merge.data.table(allDist, f25pix, all = TRUE)
  # add the growth_curve_id to help the match??
  f25dist <- f25pix[sim$allPixDT, on = 'pixelIndex', nomatch = 0][,.(pixelIndex,age2015, age2020, yrf25,
                                                                     eventf25, growth_curve_id)]
  f25dist[, targetPix := 0L]
  f25dist[pixelIndex %in% toMatch5$pixelIndex]$targetPix <- 1
  # are there any pixels with targetPix == 1 and yrf9 !is.na()?
  f25dist[targetPix >0 & !is.na(yrf25)]

  newF25dist <- f25dist[targetPix >0 & !is.na(yrf25), .(pixelIndex, yrf25, eventf25, growth_curve_id)]
  ## assign age at time of dist

  ## calculate age1985

  rndF25BygcID <- newF25dist[,.N, by = "growth_curve_id"]
  #sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex, .(pixelIndex, growth_curve_id)][,.N, by = "growth_curve_id"]

  f25agesDistDT <- sim$allPixDT[growth_curve_id %in% rndF25BygcID$growth_curve_id & ages>80,.(growth_curve_id, ages)]

  newAgeF25atDist <- f25agesDistDT[rndF25BygcID, on = "growth_curve_id", nomatch = 0]

  #newPixDT <- newPixDT[, .SD[sample(.N, N)], by = "growth_curve_id"] # take a sample of each growth_curve_id, length N
  #newAgedistPix <- newAgedistPix1[, lapply(.SD, sample(ages, size = N[1])), by = growth_curve_id]

  newAgeF25atDistPix <- newAgeF25atDist[, .SD[sample(ages,size = N[1])], by = "growth_curve_id"]

  # This is just a test -- does each growth curve id have the same # rows as N says
  newAgeF25atDistPix[, .N == N[1],by = "growth_curve_id"]

  # merge with newF9dist
  setorder(newF25dist, growth_curve_id)
  setorder(newAgeF25atDistPix, growth_curve_id)
  set(newF25dist, NULL, "newAgeF25atDist", newAgeF25atDistPix$ages)
  newF25dist[, ageF251985 := newAgeF25atDist - (yrf25 - 1985)]

  ## add the "new" dist to allDist but make sure there are no NAs
  table(newF25dist$yrf25, useNA = "ifany") # none
  table(newF25dist$eventf25, useNA = "ifany")
  # 1    2  NaN
  # 88 4080   49
  # need to replace those NAs
  distVec <- sample(1:2,length(which(is.na(newF25dist$eventf25))), replace = TRUE)
  newF25dist[is.na(eventf25)]$eventf25 <- distVec
  f25AddDist <- newF25dist[,.(pixelIndex, yrf25, eventf25)]
  setnames(f25AddDist,c("yrf25", "eventf25"), c("year", "events"))

  dim(allDist)
  allDist <- rbind(allDist, f25AddDist)

  ## add this age column to the ageNA1985 DT
  f25ages <- newF25dist[,.(pixelIndex, ageF251985)]
  ageNA1985 <- merge.data.table(ageNA1985, f25ages, on = "pixelIndex", all.x = TRUE)

  lastMatch <- toMatch5[!(pixelIndex %in% newF25dist$pixelIndex)]

  ##these will get a random age in1985 with no disturbances
  rndLastBygcID <- sim$allPixDT[pixelIndex %in% lastMatch$pixelIndex, .(pixelIndex, growth_curve_id)][,.N, by = "growth_curve_id"]

  lastAgeDistDT <- sim$allPixDT[growth_curve_id %in% rndLastBygcID$growth_curve_id & ages>80,.(growth_curve_id, ages)]

  newAgedistLast <- lastAgeDistDT[rndLastBygcID, on = "growth_curve_id", nomatch = 0]

  #newPixDT <- newPixDT[, .SD[sample(.N, N)], by = "growth_curve_id"] # take a sample of each growth_curve_id, length N
  #newAgedistPix <- newAgedistPix1[, lapply(.SD, sample(ages, size = N[1])), by = growth_curve_id]

  newAgedistLastPix <- newAgedistLast[, .SD[sample(ages,size = N[1])], by = "growth_curve_id"]

  # This is just a test -- does each growth curve id have the same # rows as N says
  newAgedistLastPix[, .N == N[1],by = "growth_curve_id"]

  # re-attach a pixelIndex
  pixIndDistPixToMatch <- sim$allPixDT[pixelIndex %in% lastMatch$pixelIndex, .(pixelIndex, growth_curve_id)]
  setorder(pixIndDistPixToMatch, growth_curve_id)
  setorder(newAgedistLastPix, growth_curve_id)
  set(pixIndDistPixToMatch, NULL, "lastAge1985", newAgedistLastPix$ages)
  pixIndDistPixToMatch[, growth_curve_id := NULL]
  ageNA1985 <- merge.data.table(ageNA1985, pixIndDistPixToMatch, on = "pixelIndex", all.x = TRUE)

  ### NOW, make the one column
  ageNA1985[, age1985 := lastAge1985]
  #table(ageNA1985[noDist == 1 & substract>0,]$age1985, useNA = "ifany")
  ageNA1985[noDist == 1 & substract >= 0,]$age1985 <- ageNA1985[noDist == 1 & substract >= 0,]$PixSubtract1985
  #table(ageNA1985[noDist == 0 & substract>0,]$age1985, useNA = "ifany")
  ageNA1985[firePix1985age>0]$age1985 <- ageNA1985[firePix1985age>0]$firePix1985age
  ageNA1985[cutPix1985age>0]$age1985 <- ageNA1985[cutPix1985age>0]$cutPix1985age
  #table(ageNA1985[noDist == 0 & substract<=0,]$age1985, useNA = "ifany")
  ageNA1985[noDist == 0 & substract<=0,]$age1985 <- ageNA1985[noDist == 0 & substract<=0,]$distPixNegAge1985
###here
  ageNA1985[!is.na(ageF91985),]$age1985 <- ageNA1985[!is.na(ageF91985),]$ageF91985
  ageNA1985[!is.na(ageF251985),]$age1985 <- ageNA1985[!is.na(ageF251985),]$ageF251985

  #table(ageNA1985$age1985, useNA = "ifany")
  age1985filled <- ageNA1985[,.(pixelIndex, age1985)]

  # add the values to the ageRaster2015
  ageRaster2015[age1985filled$pixelIndex] <- age1985filled$age1985

  # replace the age column in sim$allPixDT with the age1985 filled values
  setorder(sim$allPixDT, pixelIndex)
  # get rid of the ages (==age2020)
  sim$allPixDT[, ages := NULL]
  # add the age1985filled
  sim$allPixDT[, ages := ageRaster2015[]]

  return(invisible(sim))
}
