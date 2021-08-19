#This is an SEIQ model for simulating an outbreak of acute-symptom inducing African Swine Fever Virus (ASFV) in a Brazilian swine-production system.

#model terms:

#edgeSet = Time-ordered edgeSet derived from TERGM-initialized networks. The data frame should have 4 columns: onset, terminus, tail, and head (the same columns required for creating a temporally-dynamic network object using the statnet package). Note that the time unit is assumed to be days.
#nodeSet = Data frame containing columns related to node characteristics utilized for TERGM-creation. Note that the function assumes that a "node_id" column exists in the nodeSet.
#distMatrix = Square matrix containing inter-node distances (in km) for each node pair. Note that distances must be in kilometers! Additionally, nodes must be reported as 0-km from themselves. Aside from these, there are no other requirements for reported distances (i.e., distances can be: Euclidean distances, great-circle distances, smallest-road distances, etc.)
#virulence = Character string taking the values "moderate" or "high" that dictates the strain virulence (ie., potential time (in days) to mortality). If "moderate," mortality results in 11-15 days. If "high" mortality results in 6-9 days. Defaults to "moderate." Note that time to mortality is taken from the FAO Production and Health Manual: African Swine Fever: Detection and Diagnosis (Beltran-Alrcudo et al. 2017).
#systemSet = Matrix, data frame, or ist of matrices/data frames with columns representing production systems (i.e., companies) and a row for each node in the distMatrix. Note: if input is a list of matrices/data frames, each list entry must be representative of the systems at a given timestep (e.g., list entry 1 is for timestep 1). Note: Values within these matrices must be binary with "1" indicating system presence and "0" indicating system absence. Additionally, if nodes are associated with unknown systems, there must be an associated "UNK" column in the input.
#TERGM = TERGM (or ERGM) object to be used for network initialization.
#transProbE = Probability of transmission of an exposed individual from an infectious node to a susceptible one given no probability scaling. Note that this can be thought of as the probability of a "susceptible" node transitioning to the "exposed" state following contact with another infectious node.
#transProbI = Probability of transmission of an infectious individual from an infectious node to a susceptible one given no probability scaling. Note that this can be thought of as the probability of a "susceptible" node transitioning to the "infectious" state following contact with another "infectious" node.
#transProbScaling = Daily increase in ProbTransmission value when nodes are in the infectious state (I). This represents intra-node transmission.
#quarantineContinue = Logical. If True, quarantine will be mandated for I-state notes once the infection has been discovered (i.e., through sampling or first mortality event), and these nodes will be removed from the network and excluded from later simulation. If False, the simulation will end after the first infection is discovered. (NOTE: quarantineContine == FALSE is not currently supported, but we plan to allow this feature in future iterations of this model)
#quarantineAfterFirstDeath = Logical. If True, infectious periods (i.e., time before quarantine) for nodes will be drawn from a beta-Pert distribution. If FALSE, infectious periods are based on within-node SEIR ODE models describing the time it takes for mortalityPercentage (e.g., 0.1) percent of animals within the node to die. Note that if FALSE, the nodeSet input must contain a "numAnimals" column detailing the number of animals within each node. Defaults to TRUE.
#mortalityPercentage = Numerical. Describes the percentage of animals within nodes that must die before a quarantine is enacted if quarantineAfterFirstDeath == FALSE.
#quarProc = Integer taking the value 1, 2, 3, 4, 5, 6, or 7 that descibes the quarantine procedure to be used. If 1, nodes within quarRad distance of a quarantined nodes will also be quarantined. If 2, nodes within the same production company at the timestep of interest will be quarantined. Note that, in this case, a systemSet must be input into the model. If 3, both processes 1 and 2 take place. If 4, nodes are quarantined based on observed edges between themselves and other quarantined nodes during a defined time period (i.e., quarTimeFrame) prior to initial quarantines. If 5, both procedures 1 and 4 take place. If 6, both procedures 2 and 4 take place. If 7, procedures 1, 2, and 3 take place.
#quarRad = Radial distance (in km) of imposed quarantine around infected nodes after infection is discovered.
#quarTimeFrame = The number of temporal units (e.g., days) before the first identified infection (i.e., infection to be quarantined) during which we will mandate that edges with the quarantined node mandate additional quarantines. This is only relevant if quarProc == 4 or 5.
#originNode = The node that first transitions to state I.
#originTime  = The timestep (day) at which OriginNode transitions to state I. Note that this assumes all ties on originTime have already occurred (i.e., potential transmission on edges begin on day originTime + 1).
#terminusTime = The timestep at which the simulation ends. If NA, the simulation will proceed until no nodes are in the E or I state (NOTE: terminusTime == NA is not currently supported, but we plan to allow this feature in future iterations of this model). In this case, a new network will be initialized once all timesteps represented in the edgeSet have been simulated.

####The following terms were part of the original model, but were removed after meeting with Gustavo on 09/11/2019, as he said that random blood sampling does not take place.
#sampProb = Daily probability that farms will take blood samples from animals. 
#catchProb = Probability that at least one infected individual is sampled, if nodes are in the E or I state. Scales with transProbScale. If an infected individual is sampled, a quarantine for the affected node goes into effect DaysToID days later. 
#daysToID = A single number or numerical vector of length ≥ 2, describing the number or range of days, respectively, after which ASF will can be identified from blood tests. Note that ASF can be identified from blood tests 8-21 days post infection. 

#Additional notes:

#Note that incubation-period length ranges from 5-15 days. However, because nodes may transfer infected individuals at various incubation stages, however, the effective incubation period is equivalent to [(5-TI):15], where TI is a number between 1 and 5 representing the length of time (days) the infector has spent in the I state, and maxes out at 5 (the minimum incubation time). Similar to the incubation period, if the probability of sending infecious individuals (not just exposed one) is > 0, the same issue occurs.  
#Note that mortality of an individual occurs 6-to-13 days after the node transitions to I. When mortality occurs, a quarantine goes into effect the following day.
#Note that setting quarantineContinue == FALSE and terminusTime == NA are not currently supported. However, we intend to add this functionality to the model in the future.


#Assumptions and Limitations:

#1.) Model assumes that the modeled serotype induces acute symptoms in swine (i.e., the incubation and mortality are parameterized based on acute symptomology), not peracute or chronic symptomology. Theoretically, though, we could adjust these parameters if need be.
#2.) Model assumes ≤ 1 edge exists between nodes at each timestep.
#3.) Model assumes homogeneous probability of transmission on edges and probability scale, regardless of node characteristics.
#4.) Similarly, this model assumes that all farms sample animals at the same rate (though, I feel this is a pretty reasonable assumption to make as sampling protocols are likely specified by companies or the government), and catch infected individuals at the same rate.
#5.) Model assumes that symptomatic period coincides with infectious period.
#6.) Model assumes that each animal spends at least 28 days at nodes where they are transfered to (i.e., infected individuals have no chance of being further transfered).
#7.) Model assumes that all exposed/infectious individuals transferred to farms will die from ASF (i.e., mortality = 100%) and farms will be alerted to ASF when they do so. 
#8.) Model assumes that once a node transitions to E/I it can no longer receive infectious individuals from other nodes (i.e., lengths of time individuals are in each state are not influenced by subsequent animal transfers).
#9.) Model assumes that if quarantineAfterFirstDeath == FALSE, if no information is known about the number of animals in a grid, or the number is 0, time-to quarantine is based on a Beta-pert distribution drawn from the empirical time-to-mortality.
#10.) Along with assumptions 6&7, model assumes that infected animals will not die from anything other than ASF. This is for simplicity, but may be inaccurate (e.g., in reality slaughterhouses will likely kill animals before they die of ASF, thus they will remain unaware of ASF spread.)

ASFSim<-function(edgeSet, nodeSet, distMatrix, systemSet, TERGM = NULL, virulence = "moderate", transProbE = 1, transProbI = 0, transProbScaling = 0, quarantineContinue = F, quarantineAfterFirstDeath = TRUE, mortalityPercentage = 0.01, quarProc = 1, quarRad = 10, quarTimeFrame = 14, originNode = 2, originTime = 0, terminusTime = 1457){ ###This is a function that simulates the spreading of ASFV through a network of swine-production facilities in Brazil #Note: sampProb=0.01, catchProb= 0.05, and daysToID = seq(8,21,1) were removed 09/11/2019.
  
  daysToQuarantine.func<- function(numAnimals, infectious_period, latent_period, mortalityPercentage, transmissionRate){ #function for determining the number of days infected nodes will remain infectious until quarantined. This function assumes that nodes will be quarantined when a mortality percentage of animals die due to the infection.
    
    #Function to compute derivatives of the differential equations.
    
    seir_model = function (current_timepoint, state_values, parameters) #model is taken from the SEIR model example given by https://rpubs.com/srijana/110753
    {
      # create state variables (local variables)
      S = state_values [1]        # susceptible
      E = state_values [2]        # exposed
      I = state_values [3]        # infectious
      R = state_values [4]        # recovered (or in our case, dead)
      
      with ( 
        as.list (parameters),     # variable names within parameters can be used 
        {
          # compute derivatives
          dS = (-beta * S * I)
          dE = (beta * S * I) - (delta * E)
          dI = (delta * E) - (gamma * I)
          dR = (gamma * I)
          
          # combine results
          results = c (dS, dE, dI, dR)
          list (results)
        }
      )
    }
    
    #Compute values of beta (transmission rate) and gamma (mortality rate).
    
    beta_value = transmissionRate #pred-determined transmission rate
    gamma_value = 1 / infectious_period
    delta_value = 1 / latent_period
    
    #Disease dynamics parameters.
    
    parameter_list = c (beta = beta_value, gamma = gamma_value, delta = delta_value)
    
    #Initial values for sub-populations.
    
    W = numAnimals - 1  # susceptible hosts
    X = 1           # infectious hosts
    Y = 0           # recovered hosts
    Z = 0           # exposed hosts
    
    #Compute total population.
    
    N = numAnimals
    
    #Initial state values for the differential equations.
    
    initial_values = c (S = W/N, E = X/N, I = Y/N, R = Z/N)
    
    #Output timepoints.
    
    timepoints = seq (0, 365, by=1) #it's highly unlikely that there will be living animals after 365 days
    
    #Simulate the SEIR epidemic.
    
    output = deSolve::lsoda (initial_values, timepoints, seir_model, parameter_list)
    time = unname(output[min(which(output[,match("R", colnames(output))] >= mortalityPercentage)),1]) #pulls the time point when mortality reaches the requirement for quarantine.
    
    return(time)
    
  }
  
  #First we set up acute disease parameters (i.e., define how long individuals spend in each state.)
  
  transmissionRate <- 1.679 # the daily within-node (i.e., between individual pigs within nodes) ASF transmission rate (i.e., beta) is 1.679. This is the weighted median of estimated betas presented by Guinat et al. 2018 (title: Inferring within‐herd transmission parameters for African swine fever virus using mortality data from outbreaks in the Russian Federation), and used by Ferdousi et al. 2019 (title: Generation of swine movement network and analysis of efficient mitigation strategies for African swine fever virus) in their ASF tranmission model.
  
  E_duration.vec <- seq(4, 7, 1) #the incubation period is most-likely to be 4-7 days (Beltran-Alrcudo et al. 2017 -- FAO Production and Health Manual: African Swine Fever: Detection and Diagnosis), but will have the ability to be < 4.
  
  if(virulence == "moderate"){ #the most-likely time to mortality varies based on strain virulence.
    I_duration.vec <- seq(11,15,1) #the time to mortality is 11-15 days following onset of disease.
  }
  if(virulence == "high"){ #the time to mortality varies based on strain virulence.
    I_duration.vec <- seq(6,9,1) #the time to mortality is 6-9 days following onset of disease.
  }  
  
  #Now we set up our record-keeping objects.
  
  parameterFrame<-data.frame(virulence = virulence, transProbE = transProbE, transProbI = transProbI, transProbScaling = transProbScaling, quarantineContinue = quarantineContinue, quarProc = quarProc, quarRad = quarRad, quarTimeFrame = quarTimeFrame, quarantineAfterFirstDeath = quarantineAfterFirstDeath, mortalityPercentage = mortalityPercentage, originNode = originNode, originTime = originTime, terminusTime = terminusTime) #create a data frame recording all specified parameters to be included in output for reference purposes.
  
  node.vec <- unique(nodeSet[,match("node_id", colnames(nodeSet))]) # a vector of node ids. Note that the first column of nodeSet must reflect node IDs.
  
  epiStateTracking.list <-list(NULL) #set up a list to track the state variables over time. This will be a model output.
  listTrack <- 1 #An object identifying the number of objects included in epiStateTracking.list. This number will update at each timestep.
  
  state.matrix <- matrix(ncol = 5, nrow = length(node.vec)) #create a matrix containing transition states at each timestep
  colnames(state.matrix)<-c("node_id", "state", "transitionStep", "timeToNextState", "timeSinceTransition") #node_id is self explanatory, state is the state at the current timestep, transitionStep is the timestep at which nodes transitioned to the current state, timeToNextState is the number of days following the transitionStep required to transition to the next state, timeSinceTransition is the number of days that have passed since the state was changed 
  infection_matrix <- matrix(ncol = 9, nrow = length(node.vec)) #create a matrix to track if nodes were infected during simulation and how many times they were infected. This will be a model output.
  colnames(infection_matrix)<-c("node_id", "infected", "infectedBy", "infectedTime", "infectionDist", "quarantinedBy", "quarantinedTime", "originNode", "originDist") #node_id is self explanatory, infected is binary with "0" indicating never transitioning to E or I states and "1" indicating they transitioned at least once, "infectedBy" is a categorical variable indicated what infectious node caused the state transition, "infectionDist" is the numerical distance between the infected and infector nodes, "quarantinedBy" is a categorical variable describing the nodeID of the individual whose quarantine radius resulted in a node being removed from the simulation, "originNode" is binary indicating if nodes were chosen as the originNode in this simulation, and "originNode" is the distance between infected nodes and the origin node: used to evaluate the spread of epidemic.
  time.matrix <- matrix(ncol = 5, nrow = (1 + (terminusTime - originTime))) #create a matrix for describing the infected nodes and max spread distance (i.e., distance from the origin node) at each successive time point 
  colnames(time.matrix)<- c("timeStep", "transmissionFrom", "transmissionTo", "infectedList", "maxSpreadDistance")
  
  state.matrix[,1] <- node.vec #denote the node IDs
  state.matrix[,2] <- "S"  #all nodes begin as susceptible except for the origin node. It's easier to call them all "S" initially and immediately redefine the state of the originNode
  state.matrix[which(node.vec == originNode),2] <- "I" #redefine the state of the originNode as "infectious."
  state.matrix[,3] <- originTime #All of these states began at the originTime timestep.
  
  if(quarantineAfterFirstDeath == TRUE){ #if this parameter is true, the time before quarantine (same as the infectious period), will be obtained by sampling from a beta-pert distribution
    state.matrix[which(node.vec == originNode),4] <-round(mc2d::rpert(n = 1, min = 0, mode = mean(I_duration.vec), max = max(I_duration.vec), shape = 4)) #determine how long node will remain infectious before becoming quarantined by sampling from a beta-PERT distribution (Vose 2008). Note: value is rounded because it needs to be an integer. Additional note: infectious time is most-likely to be the mean of the infectious-time range described by Beltran-Alrcudo et al. 2017 -- FAO Production and Health Manual: African Swine Fever: Detection and Diagnosis, but will have the ability to be < min(I_duration.vec).
  }else{ #if the parameter is FALSE, the time before quarantined is determined via an intra-node mechanistic compartmental transmission model. The time to quarantine (synonomous with the infectious period) is taken from the time required for a fixed percentage of the population within the node (mortalityPercentage) to die from the disease.
    if(is.na(nodeSet[which(node.vec == originNode), match("numAnimals", colnames(nodeSet))]) == FALSE && nodeSet[which(node.vec == originNode), match("numAnimals", colnames(nodeSet))] >= 1){ #if the capacity is NA or 0, the model treats the situation as if quarantineAfterFirstDeath == TRUE.
      state.matrix[which(node.vec == originNode),4] <-daysToQuarantine.func(numAnimals = nodeSet[which(node.vec == originNode), match("numAnimals", colnames(nodeSet))], infectious_period = mean(I_duration.vec), latent_period = mean(E_duration.vec), mortalityPercentage, transmissionRate)
    }else{ #if the numAnimals is NA or 0
      state.matrix[which(node.vec == originNode),4] <-round(mc2d::rpert(n = 1, min = 0, mode = mean(I_duration.vec), max = max(I_duration.vec), shape = 4)) #determine how long node will remain infectious before becoming quarantined by sampling from a beta-PERT distribution (Vose 2008). Note: value is rounded because it needs to be an integer. Additional note: infectious time is most-likely to be the mean of the infectious-time range described by Beltran-Alrcudo et al. 2017 -- FAO Production and Health Manual: African Swine Fever: Detection and Diagnosis, but will have the ability to be < min(I_duration.vec).
    }
  }
  state.matrix[,5] <- 0 #all nodes transitioned at this timestep.
  
  infection_matrix[,1] <- node.vec #denote the node IDs
  infection_matrix[,2] <- 0 #all nodes begin as susceptible except for the origin node. It's easier to call them all 0 initially and immediately redefine the state of the originNode.
  infection_matrix[which(node.vec == originNode),2] <- 1 #redefine the state of the originNode
  infection_matrix[,3] <- NA #all nodes begin as susceptible except for the origin node, which was not infected by anyone.
  infection_matrix[,4] <- NA #all nodes begin as susceptible except for the origin node
  infection_matrix[which(node.vec == originNode),4] <- originTime # the originNode was infected at the first timestep
  infection_matrix[,5] <- NA #all nodes begin as susceptible except for the origin node, which was not infected by anyone.
  infection_matrix[,6] <- NA #no nodes have been quarantined at the first timestep
  infection_matrix[,7] <- NA #no nodes have been quarantined at the first timestep
  infection_matrix[,8] <- 0 #most nodes are not origin nodes. It's easier to call them all 0 initially and immediately redefine the state of the originNode.
  infection_matrix[which(node.vec == originNode),8] <- 1 #redefine the state of the originNode
  infection_matrix[,9] <- NA #most nodes are not infected yet. It's easier to call them all NA initially and immediately redefine the state of the originNode.
  infection_matrix[which(node.vec == originNode),9] <- 0 #redefine the state of the originNode. It is 0 km away from itself.
  
  time.matrix[,1] <- seq(originTime, terminusTime, by = 1) #create the sequence of time steps
  time.matrix[1,3] <- originNode #only the origin node transitions at the first time step
  time.matrix[1,4] <- originNode #only the origin node transitions at the first time step
  time.matrix[1,5] <- 0 #the infection has not spread away from the origin node. 
  
  epiStateTracking.list[[listTrack]] <-state.matrix #update the continuous state tracker. #Note that the first list object in this output will always reflect states at originTime
  
  current_timeStep <- originTime #simulations start on originTime 
  sequence_timeStep <- 1 #the originTime is the first time step in the timestep sequence.
  
  stop_sim <- F #logical object that identifies whether or not the simulation should be stopped (if TRUE). 
  
  while(sum(length(which(state.matrix[,2] == "I")), length(which(state.matrix[,2] == "E"))) > 0 & stop_sim == F){ #a loop is required here because we need to consistently keep up with changes over time. Note that the loop ends when no exposed or infectious nodes exist, or stop_sim is set to TRUE.
    
    listTrack <- listTrack + 1 #update the number of list objects in the continuous track.
    current_timeStep <- current_timeStep + 1 #update the timeStep.
    sequence_timeStep <- sequence_timeStep + 1 #update the sequence timeStep
    
    #Determine if any infections occur
    transfersAtTimeStep<- droplevels(subset(edgeSet, onset == current_timeStep)) # create a data frame containing transfers occurring at the current_timestep.
    stateChangeAtTimeStep.vec<-NULL #A vector containing the node IDs of state-changing individuals at currentTimeStep
    
    if(nrow(transfersAtTimeStep) > 0){ #there is a possibility that no transfers will occur on the timestep. In this case, transfersAtTimeStep will fave 0 rows and we can skip all transmission calculations.
      
      S_set.vec <- state.matrix[which(state.matrix[,2] == "S"),1] #identify which nodes are currently in the susceptible state. 
      I_set.vec <- state.matrix[which(state.matrix[,2] == "I"),1] #identify which nodes are currently in the infectious state. 
      
      potentialTransmissionEdges<-droplevels(transfersAtTimeStep[which(transfersAtTimeStep$head%in%S_set.vec == T),]) #identify which to nodes are susceptible to receiving infectious animals on the timeStep, and subset the data set accordingly.
      
      I_Transmission.vec <-I_set.vec[which(I_set.vec%in%potentialTransmissionEdges$tail == T)] #pull the IDs of infectious nodes observed to sending animals at the timestep.
      
      if(length(I_Transmission.vec) > 0){ #if any from (i.e., tail) nodes in the potentialTransmissionEdges edge set are in the I state. There's a possibility that none will exist, and we can skip the transmission calculations.
        I_evaluationOrder <- sample(I_Transmission.vec, size = length(I_Transmission.vec), replace = FALSE) #randomly decide on evaluation order of nodes. This is important because we are interested in capturing what nodes have increased probability of transmitting ASF, and we don't want node_id order to bias this estimate.
        for(i in I_evaluationOrder){ # for-loop to evaluate transmission from each node
          I_timeInState <- as.integer(state.matrix[which(state.matrix[,1] == i),5]) #determine the length of time infectious nodes have been in the I state.
          probabilityScaling<-transProbScaling*I_timeInState #The probability of transmission increases by a fixed percentage each day a node has been infectious.
          I_edges<-droplevels(subset(potentialTransmissionEdges, tail == i)) #subset edges at timestep current_timestep originating from node i. Note: if it gets to this point, there will always be ≥1 edge in edgeSet.
          potentialInfectedNodes<- I_edges$head #identify the susceptible nodes that may change states due to the influence of node i.
          
          for(j in potentialInfectedNodes){ #for-loop evaluating transmission to each specific node.
             exposed_sample <- sample(seq(0.01,1, by = 0.01), 1, replace = F) #sample from a distribution between 0.01 and 1 to determine if nodes become exposed.
             isExposed <- ifelse(exposed_sample <= ((transProbE*probabilityScaling) + transProbE), T, F) # TRUE or FALSE if node j becomes exposed.
         
            if(isExposed == T){ # if sampling indicated node should shift to state E.
                
              stateChangeAtTimeStep.vec<- c(stateChangeAtTimeStep.vec, j) #add j to the vector of state-changing nodes at this timestep.
              state.matrix[which(state.matrix[,1] == j), 2] <- "E" #change the state of node j to "E."
              state.matrix[which(state.matrix[,1] == j), 3] <- current_timeStep #change the timestep of transition to current_timestep.
              state.matrix[which(state.matrix[,1] == j), 4] <- round(mc2d::rpert(n = 1, min = 0, mode = mean(E_duration.vec), max = max(E_duration.vec), shape = 4)) #determine how long node will remain exposed before becoming infectious by sampling from a beta-PERT distribution (Vose 2008). Note: value is rounded because it needs to be an integer. Additional note: exposed time is most-likely to be the mean of the exposure-time range described by Beltran-Alrcudo et al. 2017 -- FAO Production and Health Manual: African Swine Fever: Detection and Diagnosis, but will have the ability to be < 4.
              state.matrix[which(state.matrix[,1] == j),5] <- 0 #change the timeSinceTransition to 0
              
              infection_matrix[which(infection_matrix[,1] == j), 2] <- 1 #change the infection status of node j in the infection_matrix
              infection_matrix[which(infection_matrix[,1] == j), 3] <- i #report who infected node j (i.e., node i).
              infection_matrix[which(infection_matrix[,1] == j), 4] <- current_timeStep #report the time that j was infected.
              infection_matrix[which(infection_matrix[,1] == j), 5] <- distMatrix[i,j] #report the distance between nodes i & j. 
              infection_matrix[which(infection_matrix[,1] == j), 9] <- distMatrix[which(node.vec == originNode),j] #report the distance between node j and the origin node. 
            }
          } #j for-loop
        } #i for-loop
      } #if there were infectious from-nodes
    } #if there were any animal transfers at current_timestep
    
    #Now update time in state variables, and determine if state variables change due to time spent in specific states.
    state.matrix[which(state.matrix[,1]%in%stateChangeAtTimeStep.vec == FALSE),5]<-as.integer(state.matrix[which(state.matrix[,1]%in%stateChangeAtTimeStep.vec == FALSE),5]) + 1 #add one day to the timeSinceTransition column for nodes that did not become infected earlier in the simulation.
    newINodes<-state.matrix[which(((as.integer(current_timeStep) - as.integer(state.matrix[,3])) == as.integer(state.matrix[,4])) == TRUE & state.matrix[,2] == "E"),1] #pull the ids of nodes that need to transition from E to I.
    newQNodes<-state.matrix[which(((as.integer(current_timeStep) - as.integer(state.matrix[,3])) == as.integer(state.matrix[,4])) == TRUE & state.matrix[,2] == "I"),1] #pull the ids of nodes that need to transition from I to Q.
    
    if(length(newINodes) > 0){ #If there are any new I-state nodes transitioned from E states.
      for(h in newINodes){ #for-loop to assign new state variables.
        state.matrix[which(state.matrix[,1] == h), 2] <- "I" #change the state of node j to "E."
        state.matrix[which(state.matrix[,1] == h), 3] <- current_timeStep #change the timestep of transition to current_timestep.
        if(quarantineAfterFirstDeath == TRUE){ #if this parameter is true, the time before quarantine (same as the infectious period), will be obtained by sampling from a beta-pert distribution
          state.matrix[which(state.matrix[,1] == h), 4] <- round(mc2d::rpert(n = 1, min = 0, mode = mean(I_duration.vec), max = max(I_duration.vec), shape = 4)) #determine how long node will remain infectious before becoming quarantined by sampling from a beta-PERT distribution (Vose 2008). Note: value is rounded because it needs to be an integer. Additional note: infectious time is most-likely to be the mean of the infectious-time range described by Beltran-Alrcudo et al. 2017 -- FAO Production and Health Manual: African Swine Fever: Detection and Diagnosis, but will have the ability to be < min(I_duration.vec).
        }else{ #if the parameter is FALSE, the time before quarantined is determined via an intra-node mechanistic compartmental transmission model. The time to quarantine (synonomous with the infectious period) is taken from the time required for a fixed percentage of the population within the node (mortalityPercentage) to die from the disease.
          if(is.na(nodeSet[which(node.vec == h), match("numAnimals", colnames(nodeSet))]) == FALSE && nodeSet[which(node.vec == h), match("numAnimals", colnames(nodeSet))] >= 1){ #if the capacity is NA or 0, the model treats the situation as if quarantineAfterFirstDeath == TRUE.
            state.matrix[which(state.matrix[,1] == h),4] <-daysToQuarantine.func(numAnimals = nodeSet[which(node.vec == h), match("numAnimals", colnames(nodeSet))], infectious_period = mean(I_duration.vec), latent_period = mean(E_duration.vec), mortalityPercentage, transmissionRate)
          }else{ #if numAnimals is NA or 0
            state.matrix[which(state.matrix[,1] == h), 4] <- round(mc2d::rpert(n = 1, min = 0, mode = mean(I_duration.vec), max = max(I_duration.vec), shape = 4)) #determine how long node will remain infectious before becoming quarantined by sampling from a beta-PERT distribution (Vose 2008). Note: value is rounded because it needs to be an integer. Additional note: infectious time is most-likely to be the mean of the infectious-time range described by Beltran-Alrcudo et al. 2017 -- FAO Production and Health Manual: African Swine Fever: Detection and Diagnosis, but will have the ability to be < min(I_duration.vec).
          }
        }
        state.matrix[which(state.matrix[,1] == h), 5] <- 0 #change the timeSinceTransition to 0
      }
    }
    
    if(length(newQNodes) > 0){#If there are any new Q-state nodes transitioning from I states.
      Q_evaluationOrder <- sample(newQNodes, size = length(newQNodes), replace = FALSE) #As with randomizing the order of infectious nodes, we randomly decide on evaluation order of nodes. 
      for(k in Q_evaluationOrder){ #for-loop for identifying quarantined nodes. Note that additional nodes within quarRad will also be quarantined.
        if(quarProc == 1){ #if the quarantine procedure is based on distance from I->Q transitioning nodes
          quarantine.vec <- unname(which(as.vector(distMatrix[,as.integer(k)]) <= quarRad)) #create a vector describing what nodes should be quarantined due to an infection at node h at this timestep. Note that all nodes within a quarRad distance will be quarantined. Additional note: the as.vector command is neccessary to remove any "units" associated with the matrix.
        }
        if(quarProc == 2){ #if the quarantine procedure is based on shared production system with I->Q transitioning nodes
          if(is.list(systemSet) == TRUE && is.data.frame(systemSet) == TRUE){ #if the systemSet is a non-data frame list, the function assumes there is one entry for each time step
            systemSet.timeStep<-systemSet[[current_timeStep]] #pull the systemSet for the final timeStep
          }else{ #if the function is not a list, or is a data frame, there is no need to pull anything else
            systemSet.timeStep <- systemSet
          }
          quarantine.vec<-NULL #create an empty vector to be filled below
          node_systems <- unlist(systemSet.timeStep[as.integer(k),]) #pull the system information for this node.
          system.IDs <- names(node_systems)[which(node_systems == 1)] #identify what specific production systems are associated with the node
          
          if(length(which(system.IDs == "UNK")) > 0){ #if there are unknown ("UNK") systems, we must remove that classification here. This is because quarantining all nodes with unknown systems is not reasonable.
            system.IDs <- system.IDs[-which(system.IDs == "UNK")] 
          }
          
          if(length(system.IDs) > 0){ #It's possible that removing unknowns will remove all observations. In that case (i.e., if length(system.IDs) == 0), the only node quarantine will be the k node.
            for(l in 1:length(system.IDs)){ #loop through production systems
              quarantine.vec<-c(quarantine.vec, which(systemSet[,match(system.IDs[l], colnames(systemSet))] == 1)) #bind the identified nodes to quarantine.vec
            } 
            quarantine.vec <- unique(quarantine.vec) #reduce this vector to unique values.
          }else{ 
            quarantine.vec <- k #only quarantine the k node
          } 
        }
        if(quarProc == 3){ #if the quarantine procedure is based on shared production system with I->Q transitioning nodes AND distance from I->Q transitioning nodes
          if(is.list(systemSet) == TRUE && is.data.frame(systemSet) == TRUE){ #if the systemSet is a non-data frame list, the function assumes there is one entry for each time step
            systemSet.timeStep<-systemSet[[current_timeStep]] #pull the systemSet for the final timeStep
          }else{ #if the function is not a list, or is a data frame, there is no need to pull anything else
            systemSet.timeStep <- systemSet
          }          
          quarantine.vec <- unname(which(as.vector(distMatrix[,as.integer(k)]) <= quarRad)) #create a vector describing what nodes should be quarantined due to an infection at node h at this timestep. Note that all nodes within a quarRad distance will be quarantined. Additional note: the as.vector command is neccessary to remove any "units" associated with the matrix.
          node_systems <- unlist(systemSet.timeStep[as.integer(k),]) #pull the system information for this node.
          system.IDs <- names(node_systems)[which(node_systems == 1)] #identify what specific production systems are associated with the node
          
          if(length(which(system.IDs == "UNK")) > 0){ #if there are unknown ("UNK") systems, we must remove that classification here. This is because quarantining all nodes with unknown systems is not reasonable.
            system.IDs <- system.IDs[-which(system.IDs == "UNK")] 
          }
          
          if(length(system.IDs) > 0){ #It's possible that removing unknowns will remove all observations. In that case (i.e., if length(system.IDs) == 0), the only node quarantine will be the k node.
            for(l in 1:length(system.IDs)){ #loop through production systems
              quarantine.vec<-c(quarantine.vec, which(systemSet[,match(system.IDs[l], colnames(systemSet))] == 1)) #bind the identified nodes to quarantine.vec
            } 
          }
          
          quarantine.vec <- unique(quarantine.vec) #reduce this vector to unique values.
        }
        
        if(quarProc == 4){ #if the quarantine procedure is based on previous edges between nodes
          timePeriod_transfers <- droplevels(subset(edgeSet, onset >= (current_timeStep - quarTimeFrame) & onset <= current_timeStep)) # create a data frame containing transfers occurring between the current_timestep and quarTimeFrame time units prior.
          quarantineEdges <- droplevels(subset(timePeriod_transfers, head == k | tail == k)) #pull the edges involving the quarantineNode 
          quarantine.vec <- unique(c(k, quarantineEdges$tail, quarantineEdges$head)) #create a vector describing what nodes should be quarantined due to an infection at node k at this timestep. Note that all nodes with an observed edge to the quarantineNode (regardless of direction) during the described time period will be quarantined.
        }
        
        if(quarProc == 5){ #if the quarantine procedure is based on previous edges between nodes AND quarantine radius
          quarantine.vec <- unname(which(as.vector(distMatrix[,as.integer(k)]) <= quarRad)) #create a vector describing what nodes should be quarantined due to an infection at node h at this timestep. Note that all nodes within a quarRad distance will be quarantined. Additional note: the as.vector command is neccessary to remove any "units" associated with the matrix.
          timePeriod_transfers <- droplevels(subset(edgeSet, onset >= (current_timeStep - quarTimeFrame) & onset <= current_timeStep)) # create a data frame containing transfers occurring between the current_timestep and quarTimeFrame time units prior.
          quarantineEdges <- droplevels(subset(timePeriod_transfers, head == k | tail == k)) #pull the edges involving the quarantineNode 
          quarantine.vec <- unique(c(quarantine.vec, quarantineEdges$tail, quarantineEdges$head)) #create a vector describing what nodes should be quarantined due to an infection at node k at this timestep. Note that all nodes with an observed edge to the quarantineNode (regardless of direction) during the described time period will be quarantined.
        }
        
        if(quarProc == 6){ #if the quarantine procedure is based on previous edges between nodes AND shared production system
          if(is.list(systemSet) == TRUE && is.data.frame(systemSet) == TRUE){ #if the systemSet is a non-data frame list, the function assumes there is one entry for each time step
            systemSet.timeStep<-systemSet[[current_timeStep]] #pull the systemSet for the final timeStep
          }else{ #if the function is not a list, or is a data frame, there is no need to pull anything else
            systemSet.timeStep <- systemSet
          }
          quarantine.vec<-NULL #create an empty vector to be filled below
          node_systems <- unlist(systemSet.timeStep[as.integer(k),]) #pull the system information for this node.
          system.IDs <- names(node_systems)[which(node_systems == 1)] #identify what specific production systems are associated with the node
          
          if(length(which(system.IDs == "UNK")) > 0){ #if there are unknown ("UNK") systems, we must remove that classification here. This is because quarantining all nodes with unknown systems is not reasonable.
            system.IDs <- system.IDs[-which(system.IDs == "UNK")] 
          }
          
          if(length(system.IDs) > 0){ #It's possible that removing unknowns will remove all observations. In that case (i.e., if length(system.IDs) == 0), the only node quarantine will be the k node.
            for(l in 1:length(system.IDs)){ #loop through production systems
              quarantine.vec<-c(quarantine.vec, which(systemSet[,match(system.IDs[l], colnames(systemSet))] == 1)) #bind the identified nodes to quarantine.vec
            } 
          }else{ 
            quarantine.vec <- k #only quarantine the k node
          } 
          timePeriod_transfers <- droplevels(subset(edgeSet, onset >= (current_timeStep - quarTimeFrame) & onset <= current_timeStep)) # create a data frame containing transfers occurring between the current_timestep and quarTimeFrame time units prior.
          quarantineEdges <- droplevels(subset(timePeriod_transfers, head == k | tail == k)) #pull the edges involving the quarantineNode 
          quarantine.vec <- unique(c(quarantine.vec, quarantineEdges$tail, quarantineEdges$head)) #create a vector describing what nodes should be quarantined due to an infection at node k at this timestep. Note that all nodes with an observed edge to the quarantineNode (regardless of direction) during the described time period will be quarantined.
        }
        
        if(quarProc == 7){ #if the quarantine procedure is based on previous edges between nodes AND quarantine radius AND shared production system
          if(is.list(systemSet) == TRUE && is.data.frame(systemSet) == TRUE){ #if the systemSet is a non-data frame list, the function assumes there is one entry for each time step
            systemSet.timeStep<-systemSet[[current_timeStep]] #pull the systemSet for the final timeStep
          }else{ #if the function is not a list, or is a data frame, there is no need to pull anything else
            systemSet.timeStep <- systemSet
          }          
          quarantine.vec <- unname(which(as.vector(distMatrix[,as.integer(k)]) <= quarRad)) #create a vector describing what nodes should be quarantined due to an infection at node h at this timestep. Note that all nodes within a quarRad distance will be quarantined. Additional note: the as.vector command is neccessary to remove any "units" associated with the matrix.
          node_systems <- unlist(systemSet.timeStep[as.integer(k),]) #pull the system information for this node.
          system.IDs <- names(node_systems)[which(node_systems == 1)] #identify what specific production systems are associated with the node
          
          if(length(which(system.IDs == "UNK")) > 0){ #if there are unknown ("UNK") systems, we must remove that classification here. This is because quarantining all nodes with unknown systems is not reasonable.
            system.IDs <- system.IDs[-which(system.IDs == "UNK")] 
          }
          
          if(length(system.IDs) > 0){ #It's possible that removing unknowns will remove all observations. In that case (i.e., if length(system.IDs) == 0), the only node quarantine will be the k node.
            for(l in 1:length(system.IDs)){ #loop through production systems
              quarantine.vec<-c(quarantine.vec, which(systemSet[,match(system.IDs[l], colnames(systemSet))] == 1)) #bind the identified nodes to quarantine.vec
            } 
          }
            
          timePeriod_transfers <- droplevels(subset(edgeSet, onset >= (current_timeStep - quarTimeFrame) & onset <= current_timeStep)) # create a data frame containing transfers occurring between the current_timestep and quarTimeFrame time units prior.
          quarantineEdges <- droplevels(subset(timePeriod_transfers, head == k | tail == k)) #pull the edges involving the quarantineNode 
          quarantine.vec <- unique(c(quarantine.vec, quarantineEdges$tail, quarantineEdges$head)) #create a vector describing what nodes should be quarantined due to an infection at node k at this timestep. Note that all nodes with an observed edge to the quarantineNode (regardless of direction) during the described time period will be quarantined.
        }
        
        state.matrix[which(state.matrix[,1]%in%quarantine.vec == T), 2] <- "Q" #change the state to "Q."
        state.matrix[which(state.matrix[,1]%in%quarantine.vec == T), 3] <- current_timeStep #change the transitionStep to the current timestep.
        state.matrix[which(state.matrix[,1]%in%quarantine.vec == T), 4] <- NA #change the timeToNextTransition to NA.
        state.matrix[which(state.matrix[,1]%in%quarantine.vec == T), 5] <- 0 #change the timeSinceTransition to 0.
        infection_matrix[which(infection_matrix[,1]%in%quarantine.vec == T), 6] <- k #identify the node that got nodes quarantined.
        infection_matrix[which(infection_matrix[,1]%in%quarantine.vec == T), 7] <- current_timeStep #identify the timestep when nodes were quarantined.
      }
      if(quarantineContinue == F){ #if quarantineContinue is set to false, then as soon as a node enters the Q state we tell the simulation to stop.
        stop_sim <- T #tell the while loop to stop.
      }else{#if quarantine == T
        quarantinedNodes<-state.matrix[which(state.matrix[,2] == "Q"),1] #pull the IDs of all qarantined nodes
        nodeSetforTERGM<-droplevels(nodeSet[which(nodeSet[,1]%in%quarantinedNodes == F),]) #update the nodeSet so that we may reconstitute a network from a TERGM without quarantined nodes.
        #####
        #####
        #####
        #This code section will be updated in future iterations of the model. Until this is added, setting quarantineContinue == TRUE will not work properly
        #####
        #####
        #####
      }
    }
    
    #Finally, we determine if the model fas reached terminusTime and update 
    #Note that if terminusTime is NA, then the loop will run until the last timestep reported to have an edge.
    
    if(is.na(terminusTime) == FALSE){ #if terminusTime is NA, a new network may need to be initialized.
      if(terminusTime == current_timeStep){ #if the current_timeStep reaches terminusTime, the simulation stops
        stop_sim <- T #tell the simulation loop to stop
      }
    }else{ #is terminusTime set to NA?
      #####
      #####
      #####
      #This code section will be updated in future iterations of the model. Until this is added, setting terminusTime == NA will not work properly
      #####
      #####
      #####
    }
    epiStateTracking.list[[listTrack]] <-state.matrix #update the continuous state tracker.
    
    #update the time.matrix
    spaceTrack.vec<-which(infection_matrix[,4] == current_timeStep) #identify which nodes were infected
    
    if(length(spaceTrack.vec) == 0){ #if there was no transmission during this timeStep, then all we need to do is pull the infected list and max spread distance from the previous time step.
      time.matrix[sequence_timeStep,c(4,5)] <- time.matrix[sequence_timeStep - 1,c(4,5)]
    }
    
    if(length(spaceTrack.vec) == 1){
      time.matrix[sequence_timeStep,2] <- infection_matrix[spaceTrack.vec,3]
      time.matrix[sequence_timeStep,3] <- infection_matrix[spaceTrack.vec,1]
      time.matrix[sequence_timeStep,4] <- paste(c(time.matrix[sequence_timeStep - 1,4], infection_matrix[spaceTrack.vec,1]), collapse = ", ")
      time.matrix[sequence_timeStep,5] <- max(infection_matrix[,9], na.rm = TRUE)
    }
    if(length(spaceTrack.vec) > 1){
      
      infected.vec<- infection_matrix[spaceTrack.vec,1]
      infectedBy.vec<- infection_matrix[spaceTrack.vec,3]
      
      #remove duplicate entries if they exist
      infected_noDup.vec <- infected.vec[which(duplicated(infected.vec) == FALSE)]
      infectedBy_noDup.vec <- infectedBy.vec[which(duplicated(infectedBy.vec) == FALSE)]
      
      time.matrix[sequence_timeStep,2] <- paste(infectedBy_noDup.vec, collapse = ", ")
      time.matrix[sequence_timeStep,3] <- paste(infected_noDup.vec, collapse = ", ")
      time.matrix[sequence_timeStep,4] <- paste(c(time.matrix[sequence_timeStep - 1,4], infected_noDup.vec), collapse = ", ")
      time.matrix[sequence_timeStep,5] <- max(infection_matrix[,9], na.rm = TRUE)
    }
    
    if(stop_sim == TRUE){ #if the loop is going to stop, remove all time steps that were not included in the simulation from time.matrix.
      time.matrix<-time.matrix[1:sequence_timeStep,]
    }
    
  }#while loop
  
  #compile output
  
  parameterFrame$lastTimeStep <- current_timeStep #report the timestep when the simulation ended
  parameterFrame$numInfected <-length(which(infection_matrix[,2] == 1)) #add a column to parameterFrame describing the number of individuals infected in this simulation.
  parameterFrame$epiDist <- max(as.numeric(infection_matrix[,9]), na.rm = T) #add a column reporting the maximum distance from the origin node. This indicates how far ASF was able to spread.
  parameterFrame$totalQuarantined <-length(which(is.na(infection_matrix[,6]) == F)) #add a column to parameterFrame describing the number of individuals quarantined in this simulation.
  parameterFrame$successfulQuarantine <-length(which(is.na(infection_matrix[,6]) == F & infection_matrix[,2] == 1)) #add a column to parameterFrame describing the number of individuals quarantined that were also infected in this simulation. Note: Unless the the simulation ended because terminusTime was reached, there will always be ≥ 1 successfully-quarantined node.
  parameterFrame$completeContainment <- ifelse(parameterFrame$numInfected == parameterFrame$successfulQuarantine, "contained", "not_contained") #add a column describing whether or not the simulation resulted in a complete containment.
  model_output <- list(data.frame(infection_matrix), parameterFrame) #these are what we want to output from the model.
  names(model_output) <-c("infectionMatrix","parameterFrame")
  
  return(model_output)
}
