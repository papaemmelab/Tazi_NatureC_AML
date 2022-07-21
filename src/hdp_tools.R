#################################
### hdp_tools.R #################
#################################
#
# This file gather some useful functions to work with Nicola Roberts HDP R package (https://github.com/nicolaroberts/hdp) which implements the hierarchical Dirichlet process
#

library('hdp')
#source('src/tools.R')


### Tools to launch HDP and get results

launch <- function(data,base_dist,initial_clusters,burn,posterior_samples,chains,aa,ab){
    number_of_chains <- chains
    chain_list <- vector('list', number_of_chains)
    hdp <- initialise_my_hdp(data = data,hh=base_dist,alphaa = aa,alphab= ab)


    for (i in 1:number_of_chains) {
        seed <- i * 100
        print_and_flush(sprintf('### Experiment %d (seed = %d) ###\n', i, seed))

        # run single hdp chain
        chain_list[[i]] <- activate_and_run_hdp(hdp,
                                                initcc = initial_clusters,
                                                burnin = burn,
                                                n      = posterior_samples,
                                                space  = 20,
                                                seed   = seed)
        print_and_flush('\n')
    }

    multi_output <- hdp_multi_chain(chain_list)
    print(multi_output)
}



### Tools to get max and second max proba with corresponding components

add_first_second_predicted_component <- function(hdp_output, data) {
    # Return a dataframe giving for each patient the probability by component, as well as the assigned component and its probability (max_proba)
    # → Arguments
    #     - hdp_output: hdpSampleChain or hdpSampleMulti object
    #     - data      : original data

    # keep all DP but first (first level)
    dd_predicted <- data.frame(comp_dp_distn(hdp_output)$mean[-1,])

# change categories colnames
    colnames(dd_predicted) <- paste0('component_', 0:(ncol(dd_predicted)-1))
    components_colnames <- colnames(dd_predicted)

    # pprint various info
    print_and_flush(sprintf('Number of components: %d\n', ncol(dd_predicted) - 1))
    print_and_flush(sprintf('Number of NA rows   : %d\n', nrow(dd_predicted[rowSums(is.na(dd_predicted)) != 0,])))

    # evaluate for each row the predicted component
    dd_predicted['initial_predicted_component'] <- apply(dd_predicted, 1, function(x) { if (all(is.na(x)))
                                                                     return(NaN)
                                                                 else
                                                                     return(which.max(x)-1)
                                                                })
    dd_predicted[, 'initial_predicted_component'] <- factor(dd_predicted[, 'initial_predicted_component'])


    dd_predicted['second_predicted_component'] <- apply(dd_predicted[,components_colnames], 1, function(x) { if (all(is.na(x)))
                                                                     return(NaN)
                                                                 else
                                                                     tmp <- names(which.max(x[x!=max(x)]))
                                                                     return(as.numeric(substr(tmp,nchar(tmp),nchar(tmp))))
                                                                })
    dd_predicted[, 'second_predicted_component'] <- factor(dd_predicted[, 'second_predicted_component'])
    
    # evaluate for each row the maximum probability associated to the predicted component
    
    dd_predicted['max_proba'] <- apply(dd_predicted[,components_colnames], 1, function(x) { if (all(is.na(x)))
                                                                                  return(NaN)
                                                                              else
                                                                                 return(max(x))
                                                                            })

    dd_predicted['second_max_proba'] <- apply(dd_predicted[,components_colnames], 1, function(x) { if (all(is.na(x)))
                                                                                  return(NaN)
                                                                              else
                                                                                 return(max(x[x!=max(x)]))
                                                                            })

    return (dd_predicted)
}


prepare_distributions <- function(df) {
    
    
    
    num_cols = ncol(df)
    bin <- function(x){
        set.seed(123)
      (rbinom(1, num_cols, mean(x))+1)/num_cols
    }

    ###Normal

    normal <- function(x){
        set.seed(123)
      abs(rnorm(1,mean(x),sd(x)))
    }

    ###Poisson

#     poisson <- function(x){
#         set.seed(123)
#       (rpois(num_cols,1))/num_cols
#     }

    ###Uniform equally over all columns

    equally <- function(x){
        set.seed(123)
      1/num_cols
    }

    ###Repet 1

    repet <- function(x){
        set.seed(123)
      1
    }
    
    
    dist_list <- NULL
    dist_list[["binomial"]] <- unlist(sapply(df,bin))
    dist_list[["gaussian"]] <- unlist(sapply(df,normal))
#     dist_list["poisson"] <- as.numeric(unlist(sapply(df,poisson)))
    dist_list[["unif"]] <- unlist(sapply(df,equally))
    dist_list[["repetition"]] <- unlist(sapply(df,repet))
    
    return(dist_list)

}



initialise_my_hdp <- function(data,hh,alphaa,alphab) {
  # Initialise a custom HDP structure (uniform base DP distribution over all categories) and assign the data on it, then return this structure as a hdpState object
  # → Arguments
  #     - data: dataframe with one row per patient and one column per category
  
  cat(sprintf('Initialise HDP on a %d x %d dataframe\n', nrow(data), ncol(data)))
  
  
  # number of categories
  n <- ncol(data)
  
  
  # initialise hdpState structure (first level)
  print_and_flush('  → create HDP structure...')
  hdp <- hdp_init(ppindex = 0,           # index of the initial DP
                  cpindex = 1,           # index of alphaa and alphab for initial DP
                  hh      = hh, # parameters of the base Dirichlet distribution: we choose uniform Dirichlet over n categories (n = ncol(data))
                  alphaa  = alphaa,           # shape hyperparameters for the gamma priors over the DP concentration parameters
                  alphab  = alphab)           # rate  hyperparameters for the gamma priors over the DP concentration parameters
  print_and_flush(' done!\n')
  
  
  # add DP nodes (second level)
  print_and_flush('  → add DP node for each patient...')
  hdp <- hdp_adddp(hdp,
                numdp = nrow(data), # add one DP for every sample in that cancer type
                ppindex= 1,          # index of parental process for the new DPs
                cpindex= 1)          # index of alphaa and alphab for each DP
  print_and_flush(' done!\n')
  
  
  # assign the data from each patient to a child DP (ie each DP in the second level receives a patient)
  print_and_flush('  → assign the data to the nodes...')
  hdp <- hdp_setdata(hdp,
                     dpindex = (1:nrow(data)) + 1, # indices of the DPs to assign data to (in same order as rows of data): we asssign data to each DP in the second level (2, ..., n + 1)
                     data = data)
  print_and_flush(' done!\n')
  
  
  return (hdp)
}
test_assignment <- function(x, dd_predicted, components_sd){
    if(x['second_predicted_component']%in%dd_predicted$predicted_component) {
        return( x['max_proba'] > x['second_max_proba']+components_sd[[x['second_predicted_component']+1]])
    }
    else {return(TRUE)} 
}
hdp_my_gridsearch <- function( data, posterior_samples, initial_clusters, burnin, chains, space, base_dist, alphaa, alphab){  #+3 components

    result_table <- data.frame('chains' = integer(),
                              'inicc' = integer(),
                              'n' = integer(),
                              'burnin' = integer(),
                              'space' = integer(),
                              'base_dist' = character(), 
                              'alphaa' = integer(),
                              'alphab' = integer(),
                              'n_components' = integer(),
                              'component_0' = integer(),
                              'assignment' = numeric(),
                              'Davies_Bouldin' =numeric(),
                              'Silhouette' = numeric(),stringsAsFactors=FALSE)
    
    
    for (a in chains) {
        for (b in initial_clusters) {
            for (c in posterior_samples) {
                for (d in burnin) {
                    for (e in space) {
                        for (f in base_dist){
                            for (g in alphaa){
                                for (h in alphab){
                                    str <-""
                                    if (f[1] == binomial[1]){
                                        str <-"binomial"
                                    } else if (f[1] == gaussian[1]){
                                        str <-"gaussian"
                                    } else if(f[1] == unif[1]){
                                        str <-"uniform"
                                    } else{
                                        str <-"repetition"
                                    }
                                    cat(sprintf('### %d chains | %d initial clusters | %d posterior samples | %d burn-in | %d space ###\n', a, b, c, d, e))
                                    flush.console()

                                    number_of_chains <- a
                                    chain_list <- vector('list', number_of_chains)
                                    hdp <- initialise_my_hdp(data=data,alphaa=g,alphab=h,hh=f)
                                    for (i in 1:number_of_chains) {
                                        seed <- i * 100

                                        activated_hdp <- dp_activate(hdp,
                                            dpindex = 1:numdp(hdp), # indices of the DP to activate: we choose all DP (1, ..., n + 1)
                                            initcc  = b,           # number of starting clusters (every data item is randomly assigned to a cluster to start with)
                                            seed    = seed)

                                        hdp_output <- hdp_posterior(activated_hdp,
                                            burnin = d, # number of burn-in iterations
                                            n      = c,      # number of posterior samples to collect
                                            space  = e,  # number of iterations between collected samples
                                            # cpiter = 1,    # The number of iterations of concentration parameter sampling to perform after each iteration.
                                            seed   = seed)

                                        chain_list[[i]] <- hdp_output

                                    }

                                    multi_output <- hdp_multi_chain(chain_list)
                                    multi_output <- hdp_extract_components(multi_output)


                                    dd_predicted <- data.frame(comp_dp_distn(multi_output)$mean[-1,])

                                    # change categories colnames
                                    colnames(dd_predicted) <- paste0('component_', 0:(ncol(dd_predicted)-1))
                                    components_colnames <- colnames(dd_predicted)

                                    # evaluate for each row the predicted component
                                    dd_predicted['predicted_component'] <- apply(dd_predicted, 1, function(x) { if (all(is.na(x)))
                                                                                                     return(NaN)
                                                                                                else
                                                                                                     return(which.max(x)-1)
                                                                                                })




                                    # evaluate for each row the maximum probability associated to the predicted component
                                    dd_predicted['max_proba'] <- apply(dd_predicted[,components_colnames], 1, function(x) { if (all(is.na(x)))
                                                                                                     return(NaN)
                                                                                                else
                                                                                                     return(max(x))})


                                    # get second predicted component
                                    n_components <- ncol(dd_predicted) - 2 - 1 # without columns 'component_0', 'predicted_component' and 'max_proba'

                                    dd_predicted['second_max_proba'] <- apply(dd_predicted[,0:n_components + 1], 1, function(x) { if (all(is.na(x)))
                                                                                                     return(NaN)
                                                                                                else
                                                                                                     return(sort(x,partial=n_components)[n_components])})

                                    dd_predicted['second_predicted_component'] <- apply(dd_predicted, 1, function(x) { if (all(is.na(x)))
                                                                                                     return(NaN)
                                                                                                else
                                                                                                     return(which(x==x['second_max_proba'])[1]-1)})



                                    # get standard deviations of components
                                    components_sd <- list()
                                    for( j in as.integer(levels(factor(dd_predicted$predicted_component)))){
                                        if(! is.na(j)){
                                            components_sd[[j+1]] <- sd(dd_predicted$max_proba[!is.na(dd_predicted$predicted_component) & dd_predicted$predicted_component == j])
                                        }
                                    }


                                    # determine which patients are well assigned
                                    dd_predicted['assignment'] <- apply(dd_predicted[, c('max_proba', 'second_max_proba','second_predicted_component')], 1, function(x) { if (all(is.na(x)))
                                                                                                     return(NaN)
                                                                                                else
                                                                                                     return(test_assignment(x, dd_predicted, components_sd))})

                                    # validation indexes
                                    traj <- data.frame(data[rowSums(data) !=0,])
                                    traj <- as.matrix(apply(traj, 2, function (x){as.numeric(x)}))
                                    part <- as.integer(dd_predicted$predicted_component[!is.na(dd_predicted$predicted_component)])
                                    intIdx <- intCriteria(traj, part, c("Davies_bouldin","Silhouette"))

                                    # result
                                    result_table[nrow(result_table)+1,] <- c(a, b, c, d, e, str, g, h, numcomp(multi_output), nrow(dd_predicted[!is.na(dd_predicted$predicted_component) &dd_predicted$predicted_component==0,]), nrow(dd_predicted[!is.na(dd_predicted$assignment) & dd_predicted$assignment==TRUE,]),intIdx$davies_bouldin,intIdx$silhouette)

                                    cat('Done!\n')
                                    flush.console()
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return(result_table)

}
initialise_hdp <- function(data) {
    # Initialise a custom HDP structure (uniform base DP distribution over all categories) and assign the data on it, then return this structure as a hdpState object
    # → Arguments
    #     - data: dataframe with one row per patient and one column per category

    cat(sprintf('Initialise HDP on a %d x %d dataframe\n', nrow(data), ncol(data)))


    # number of categories
    n <- ncol(data)

    
    # initialise hdpState structure (first level)
    print_and_flush('  → create HDP structure...')
    hdp <- hdp_init(ppindex = 0,           # index of the initial DP
                    cpindex = 1,           # index of alphaa and alphab for initial DP
                    hh      = rep(1/n, n), # parameters of the base Dirichlet distribution: we choose uniform Dirichlet over n categories (n = ncol(data))
                    alphaa  = 1,           # shape hyperparameters for the gamma priors over the DP concentration parameters
                    alphab  = 1)           # rate  hyperparameters for the gamma priors over the DP concentration parameters
    print_and_flush(' done!\n')

    
    # add DP nodes (second level)
    print_and_flush('  → add DP node for each patient...')
    hdp <- hdp_adddp(hdp,
                     numdp = nrow(data), # add one DP for every sample in that cancer type
                     pp    = 1,          # index of parental process for the new DPs
                     cp    = 1)          # index of alphaa and alphab for each DP
    print_and_flush(' done!\n')

    
    # assign the data from each patient to a child DP (ie each DP in the second level receives a patient)
    print_and_flush('  → assign the data to the nodes...')
    hdp <- hdp_setdata(hdp,
                       dpindex = (1:nrow(data)) + 1, # indices of the DPs to assign data to (in same order as rows of data): we asssign data to each DP in the second level (2, ..., n + 1)
                       data = data)
    print_and_flush(' done!\n')


    return (hdp)
}


activate_and_run_hdp <- function(hdp, initcc, burnin, n, space, seed) {
    # Activate all DP of a given HDP structure and run it, return the output as a hdpSampleChain object
    # → Arguments
    #     - hdp   : hdpState object
    #     - initcc: # number of starting clusters (every data item is randomly assigned to a cluster to start with)
    #     - burnin: number of burn-in iterations
    #     - n     : number of posterior samples to collect
    #     - space : number of iterations between collected samples
    #     - seed  : random seed for reproducible result

    cat('Activate HDP nodes and run posterior sampling\n')

    # activate the DPs to be included in the posterior sampling process (all DPS, first and second level)
    print_and_flush('  → activate HDP nodes...')
    activated_hdp <- dp_activate(hdp,
                                 dpindex = 1:numdp(hdp), # indices of the DP to activate: we choose all DP (1, ..., n + 1)
                                 initcc  = initcc,           # number of starting clusters (every data item is randomly assigned to a cluster to start with)
                                 seed    = seed)
    print_and_flush(' done!\n')

    # Run a Gibbs sampler over the activated DP nodes of a Hierarchichal Dirichlet Process, each iteration re-assigns the cluster allocation of every data item
    # Run burnin iterations, and then collect n samples from the chain with space iterations between each collected sample
    print_and_flush('  → run posterior sampling...\n')
    hdp_output <- hdp_posterior(activated_hdp,
                                burnin = burnin, # number of burn-in iterations
                                n      = n,      # number of posterior samples to collect
                                space  = space,  # number of iterations between collected samples
                                # cpiter = 1,    # The number of iterations of concentration parameter sampling to perform after each iteration.
                                seed   = seed)
    #print_and_flush(' done!\n')

    return (hdp_output)    
}


plot_posterior_sampling_chain_quality <- function(hdp_output, width = 20, height = 5) {
    # Plot the three diagnostic pplots for HDP posterior sampling chain
    # → Arguments
    #     - hdp_output: hdpSampleChain object

    # diagnostic plots for HDP posterior sampling chain
    set_notebook_plot_size(width, height)
    par(mfrow = c(1,3))

    plot_lik(hdp_output, bty = 'L')
    plot_numcluster(hdp_output, bty = 'L')
    plot_data_assigned(hdp_output, bty = 'L')
}


extract_components <- function(hdp_output) {
    # Extract components from HDP output, return a hdpSampleChain object
    # → Arguments
    #     - hdp_output: hdpSampleChain or hdpSampleMulti object

    cat('Extract HDP components from posterior sampling\n')

    print_and_flush('  → extract components...')
    hdp_output <- hdp_extract_components(hdp_output)
                                         # cos.merge = 0.9, # merge components with cosine similarity above this threshold
                                         # min.sample = 1)  # components must have more than this many samples
    print_and_flush(' done!\n')

    print_and_flush(sprintf('* %d components found\n', numcomp(hdp_output)))

    return (hdp_output)                                 
}


plot_components_size <- function(hdp_output, width = 8, height = 5) {
    # Plot the extracted components
    # → Arguments
    #     - hdp_output: hdpSampleChain or hdpSampleMulti object

    set_notebook_plot_size(width, height)
    plot_comp_size(hdp_output, bty = 'L', lab = c(3, 5, 7))
}


plot_category_distribution_by_component <- function(hdp_output, cat_names) {
    # Plot the category distribution by component
    # → Arguments
    #     - hdp_output: hdpSampleChain or hdpSampleMulti object
    #     - cat_names : names of the categories

    set_notebook_plot_size(22, 5)
    plot_comp_distn(hdp_output,
                    cat_names  = cat_names,
                    col        = "skyblue3",
                    col_nonsig = "black")
}


get_category_distribution_by_component <- function(hdp_output, cat_names, add_std = TRUE) {
    # Return a dataframe giving the category mean score and std score for each component
    # → Arguments
    #     - hdp_output: hdpSampleChain or hdpSampleMulti object
    #     - cat_names : names of the categories
    #     - add_std   : if TRUE add to the result the std score

    category_distribution <- comp_categ_distn(hdp_output)$mean

    number_of_components <- nrow(category_distribution)

    colnames(category_distribution) <- cat_names
    rownames(category_distribution) <- sprintf('mean_%d', 0:(number_of_components-1))
    component_names <- rownames(category_distribution)

    if (add_std) {
        std_by_component <- comp_categ_distn(multi_output)$cred.int

        for (i in 1:length(std_by_component)) {
            
            std <- std_by_component[[i]]
            
            colnames(std) <- cat_names
            rownames(std) <- c(sprintf('std_low_%d', i), sprintf('std_high_%d', i))

            category_distribution <- rbind(category_distribution, std)
        }
    }

    return (category_distribution)
}


get_top_categories_by_component <- function(hdp_output, cat_names, top_number) {
    # Return a dataframe giving the top categories by component
    # → Arguments
    #     - hdp_output: hdpSampleChain or hdpSampleMulti object
    #     - cat_names : names of the categories
    #     - top_number: number of categories to keep

    category_distribution <- get_category_distribution_by_component(hdp_output, cat_names, add_std = FALSE)
    top_categories <- apply(category_distribution, 1, function(x) names(sort(rank(x), decreasing = TRUE))[1:top_number])

    colnames(top_categories) <- sprintf('component_%d', 0:(ncol(top_categories)-1))

    return (top_categories)
}


get_prediction_result_dataframe <- function(hdp_output, data) {
    # Return a dataframe giving for each patient the probability by component, as well as the assigned component and its probability (max_proba)
    # → Arguments
    #     - hdp_output: hdpSampleChain or hdpSampleMulti object
    #     - data      : original data

    # keep all DP but first (first level)
    dd_predicted <- data.frame(comp_dp_distn(hdp_output)$mean[-1,])

    # change categories colnames
    colnames(dd_predicted) <- paste0('component_', 0:(ncol(dd_predicted)-1))
    components_colnames <- colnames(dd_predicted)

    # pprint various info
    print_and_flush(sprintf('Number of components: %d\n', ncol(dd_predicted) - 1))
    print_and_flush(sprintf('Number of NA rows   : %d\n', nrow(dd_predicted[rowSums(is.na(dd_predicted)) != 0,])))

    # evaluate for each row the predicted component
    dd_predicted['predicted_component'] <- apply(dd_predicted, 1, function(x) { if (all(is.na(x)))
                                                                     return(NaN)
                                                                 else
                                                                     return(which.max(x)-1)
                                                                })
    dd_predicted[, 'predicted_component'] <- factor(dd_predicted[, 'predicted_component'])

    # evaluate for each row the maximum probability associated to the predicted component
    dd_predicted['max_proba'] <- apply(dd_predicted[,components_colnames], 1, function(x) { if (all(is.na(x)))
                                                                                  return(NaN)
                                                                              else
                                                                                 return(max(x))
                                                                            })

    return (dd_predicted)
}


plot_assignement_probability_by_component <- function(dd_predicted) {
    # Plot boxplot comparing the assignement probability for each component
    # → Arguments
    #     - dd_predicted: dataframe with predicted component and max probability predicted for each patient

    set_notebook_plot_size(10, 5)

    n_initial <- nrow(dd_predicted)
    dd_predicted <- na.omit(dd_predicted)

    ggplot(dd_predicted) + theme_bw() +
        geom_boxplot(aes(x = predicted_component, y = max_proba), notch = TRUE,fill="lightsteelblue") +
        ggtitle(sprintf('Assignment probability by component (removed %d NA values)', n_initial - nrow(dd_predicted)))
}


plot_assignement_probability_distribution_by_component <- function(dd_predicted, width = 20, height = 10) {
    # Plot density comparing the assignement probability, for each component
    # → Arguments
    #     - dd_predicted: dataframe with predicted component and max probability predicted for each patient
    #     - width
    #     - height

    set_notebook_plot_size(width, height)

    n_initial <- nrow(dd_predicted)
    dd_predicted <- na.omit(dd_predicted)

    # we transform (extand) distn in a long data frame by gathering the categories columns in one column
    long_dd_predicted <- melt(dd_predicted, id = c('predicted_component', 'max_proba'))

    
    ggplot(long_dd_predicted) +
        geom_density(aes(value, fill = variable), alpha = 0.4) +
        facet_wrap(~ predicted_component, ncol = 4, labeller = label_both) + theme_bw() +
        coord_cartesian(ylim = c(0, 30)) +
        ggtitle(sprintf('Assignment probability density by component (removed %d NA values)', n_initial - nrow(dd_predicted)))
}
