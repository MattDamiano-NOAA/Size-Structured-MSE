########### Movement #################
# columns represent movement to each region

# Season 1 movement matrix (Winter)
range.vany       <- c(0.0,0.0,0.0,0.0,0.0,0.4,0.6) # VA to Montauk, NY # all have gone offshore
range.nonc       <- c(0.1,0.0,0.0,0.0,0.0,0.6,0.3) # Northern NC. # residual fish from southern regions move N and offshore
range.wmga       <- c(0.1,0.1,0.0,0.0,0.0,0.5,0.3) # Wilmington NC to East. FL
range.sofl       <- c(0.1,0.1,0.1,0.0,0.0,0.4,0.0) # Southern FL
range.carr       <- c(0.0,0.0,0.0,0.0,0.6,0.1,0.3) # Caribbean (general) # fish arrive in CAR
range.dist       <- c(0.0,0.0,0.0,0.0,0.0,0.1,0.9) # Northeast distant waters # fish move south
range.sarg       <- c(0.0,0.0,0.0,0.0,0.7,0.0,0.3) # Sargasso sea # fish move into CAR
movement.matrix1 <- rbind(range.vany, range.nonc, range.wmga, range.sofl, range.carr, range.dist, range.sarg)

# Season 2 movement matrix (Spring)
range.vany       <- c(1.0,0.0,0.0,0.0,0.0,0.0,0.0) # assume any fish that come here stay
range.nonc       <- c(0.3,0.7,0.0,0.0,0.0,0.0,0.0) # Northern NC
range.wmga       <- c(0.2,0.3,0.5,0.0,0.0,0.0,0.0) # Wilmington NC to GA
range.sofl       <- c(0.0,0.2,0.3,0.5,0.0,0.0,0.0) # Southern FL # assume half stay and half move north
range.carr       <- c(0.0,0.0,0.4,0.1,0.5,0.0,0.0) # Caribbean (general) # assume CAR fish don't go to FL as much
range.dist       <- c(0.0,0.0,0.0,0.0,0.5,0.0,0.5) # Northeast distant waters # send all NED fish to CAR and SAR
range.sarg       <- c(0.0,0.2,0.3,0.0,0.5,0.0,0.0) # Sargasso sea # send all SAR fish to CAR and north
movement.matrix2 <- rbind(range.vany, range.nonc, range.wmga, range.sofl, range.carr, range.dist, range.sarg)

# Season 3 movement matrix (Summer)
range.vany       <- c(0.8,0.0,0.0,0.0,0.0,0.2,0.0) # assume most fish stay
range.nonc       <- c(0.3,0.7,0.0,0.0,0.0,0.0,0.0) # Northern NC # same as spring
range.wmga       <- c(0.2,0.3,0.5,0.0,0.0,0.0,0.0) # Wilmington NC to GA # same as spring
range.sofl       <- c(0.0,0.3,0.4,0.3,0.0,0.0,0.0) # Southern FL # assume most go north in summer
range.carr       <- c(0.2,0.3,0.3,0.1,0.1,0.0,0.0) # Caribbean (general) # assume most move north
range.dist       <- c(0.0,0.0,0.0,0.0,0.0,0.0,1.0) # Northeast distant waters # send all NED fish and SAR
range.sarg       <- c(0.1,0.2,0.4,0.2,0.1,0.0,0.0) # Sargasso sea # send all SAR fish to CAR and north
movement.matrix3 <- rbind(range.vany, range.nonc, range.wmga, range.sofl, range.carr, range.dist, range.sarg)

# Season 4 movement matrix (Fall)
range.vany       <- c(0.7,0.0,0.0,0.0,0.0,0.3,0.0) # assume most fish stay
range.nonc       <- c(0.7,0.3,0.0,0.0,0.0,0.0,0.0) # Northern NC
range.wmga       <- c(0.5,0.5,0.0,0.0,0.0,0.0,0.0) # Wilmington NC to GA
range.sofl       <- c(0.2,0.3,0.5,0.0,0.0,0.0,0.0) # Southern FL # assume most go north in summer
range.carr       <- c(0.0,0.0,0.1,0.1,0.8,0.0,0.0) # Caribbean (general) # most fish arrive and stay, few move
range.dist       <- c(0.0,0.0,0.0,0.0,0.0,0.0,1.0) # Northeast distant waters # send all NED fish and SAR
range.sarg       <- c(0.0,0.0,0.0,0.0,1.0,0.0,0.0) # Sargasso sea # send all SAR fish to CAR
movement.matrix4 <- rbind(range.vany, range.nonc, range.wmga, range.sofl, range.carr, range.dist, range.sarg)
