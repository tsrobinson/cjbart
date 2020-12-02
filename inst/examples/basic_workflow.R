set.seed(89)

subjects <- 500
rounds <- 5
obs <- subjects*rounds

fake_data <- data.frame(A = sample(c("a1","a2","a3"), obs, replace = TRUE),
                        B = sample(c("b1","b2","b3"), obs, replace = TRUE),
                        C = sample(c("c1","c2","c3"), obs, replace = TRUE),
                        D = sample(c("d1","d2","d3"), obs, replace = TRUE),
                        E = sample(c("e1","e2","e3"), obs, replace = TRUE),

                        covar1 = rep(runif(subjects, 0 ,1),
                                     each = rounds),
                        covar2 = rep(sample(c(1,0),
                                            subjects,
                                            replace = TRUE),
                                     each = rounds),

                        id1 = rep(1:subjects, each=rounds),
                        stringsAsFactors = TRUE)

fake_data$Y <- ifelse(fake_data$E == "e2",
                      rbinom(obs, 1, fake_data$covar1),
                      sample(c(0,1), obs, replace = TRUE))

cj_model <- cjbart(data = fake_data,
                   Y_var = "Y",
                   id_var = "id1")

het_effects <- cjbart::OMCE(data = fake_data,
                            model = cj_model,
                            attribs = c("A","B","C","D","E"),
                            ref_levels = c("a1","b1","c1","d1","e1"),
                            Y_var = "Y",
                            id_var = "id1",
                            cores = 1)

summary(het_effects)
plot(het_effects, covar = "covar1")
