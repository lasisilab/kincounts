# =============================================================================
# Utility Functions for Kincount Simulations and Analysis
# =============================================================================

#' Estimate Fertility and Sibling Distribution Parameters
#'
#' Fits a zero-inflated negative binomial (ZINB) model to fertility count data
#' and derives parameters for both fertility and sibling distributions.
#'
#' @param fert_count Data frame containing a 'num_children' column with fertility counts
#' @return Data frame with distribution parameters (pi, size, prob, mu, mean, var)
#'         for both fertility (ZINB) and sibling (NB) distributions
estimate_fertility_parameters <- function(fert_count) {
  # ===== Fertility Distribution (Zero-Inflated Negative Binomial) =====
  # Fit ZINB model: num_children ~ ZINB(mu, size, pi)
  fert_fit <- pscl::zeroinfl(num_children ~ 1 | 1, data = fert_count, dist = "negbin")
  fert_coef <- fert_fit$coefficients
  
  # Extract parameters from the fitted model
  fert_mu <- exp(unname(fert_coef$count))        # NB mean (on log scale, so exponentiate)
  fert_pi <- plogis(unname(fert_coef$zero))      # Zero-inflation probability (on logit scale)
  fert_size <- fert_fit$theta                     # NB dispersion parameter
  fert_prob <- fert_size / (fert_size + fert_mu)  # Alternative NB parameterization
  
  # Calculate overall mean and variance accounting for zero-inflation
  fert_mean <- (1 - fert_pi) * fert_mu
  fert_var <- (1 - fert_pi) * fert_mu + (1 - fert_pi) * fert_mu^2 / fert_size + fert_pi * (1 - fert_pi) * fert_mu^2

  # ===== Sibling Distribution (Negative Binomial) =====
  # Derived from fertility distribution using size-biased sampling
  # (individuals with more siblings are more likely to be sampled)
  sib_pi <- 0  # No zero-inflation for siblings
  sib_mu <- fert_mu * (1 + fert_size) / fert_size   # Adjusted mean
  sib_size <- fert_size + 1                          # Adjusted dispersion
  sib_prob <- sib_size / (sib_size + sib_mu)         # Alternative parameterization
  sib_mean <- sib_mu
  sib_var <- sib_mu + sib_mu^2 / sib_size            # NB variance formula            # NB variance formula

  # Combine parameters into a summary data frame
  parameters <- data.frame(distribution = c("zero-inflated negative binomial", "negative binomial"),
                           pi = c(fert_pi, sib_pi),
                           size = c(fert_size, sib_size),
                           prob = c(fert_prob, sib_prob),
                           mu = c(fert_mu, sib_mu),
                           mean = c(fert_mean, sib_mean),
                           var = c(fert_var, sib_var))
  rownames(parameters) <- c("fertility", "sibling")
  return(parameters)
}


#' Calculate Genetic Match Probability
#'
#' Computes the probability that two individuals separated by g meioses share
#' at least min_num_seg DNA segments of length m cM or longer.
#'
#' @param g Number of meioses separating the relatives (g=2 for first cousins)
#' @param m Minimum segment length in centiMorgans (default: 6 cM)
#' @param min_num_seg Minimum number of segments required for a match (default: 2)
#' @param num_chrs Number of chromosomes (default: 22 autosomes)
#' @param genome_size Total genome size in Morgans (default: 35)
#' @return Probability of sharing at least min_num_seg segments of length >= m cM
genetic_match <- function(g, m = 6, min_num_seg = 2, num_chrs = 22, genome_size = 35) {
  m = m / 100  # Convert cM to Morgans (fractions)
  
  # Probability of sharing one segment of length m between relatives separated by g meioses
  # Formula: f = exp(-2*g*m) / 2^(2*g-2)
  f = exp(-2 * g * m) / 2^(2 * g - 2)
  
  # Total number of opportunities for segments across all chromosomes
  # Approximate as: num_chrs + genome_size * 2 * g
  n_opportunities = num_chrs + genome_size * 2 * g
  
  # Probability of sharing at least min_num_seg segments
  # Using binomial model: P(X >= min_num_seg) = 1 - P(X < min_num_seg)
  pr = 1 - pbinom(min_num_seg - 1, n_opportunities, f)
  
  return(pr)
}

#' Calculate Identity Probability via Cousin Matches
#'
#' Estimates the probability of identifying an individual by finding at least
#' one genetic match among their cousins in a database.
#'
#' @param relative_counts Data frame with columns: num_first_cousins, num_second_cousins, num_third_cousins
#' @param pop_size Total population size
#' @param db_size Number of individuals in the genetic database
#' @return Named vector with average probabilities for each cousin degree
identity_probability <- function(relative_counts, pop_size, db_size) {
  # Extract cousin counts from the data
  num_first_cousins <- relative_counts$num_first_cousins
  num_second_cousins <- relative_counts$num_second_cousins
  num_third_cousins <- relative_counts$num_third_cousins
  num_cousins <- cbind(num_first_cousins, num_second_cousins, num_third_cousins)
  
  # Calculate number of meioses for each cousin degree
  # First cousins: g=2, Second cousins: g=3, Third cousins: g=4
  gs <- c(1, 2, 3) + 1
  
  # Probability of genetic match for each cousin degree
  pr_genetic_match <- sapply(gs, genetic_match)
  
  # Probability that a specific cousin is in the database AND genetically detectable
  pr_detection <- pr_genetic_match * (db_size / pop_size)
  
  # For each simulation and cousin degree, calculate probability of at least one match
  # P(at least 1 match) = 1 - P(no matches) = 1 - (1 - pr_detection)^n
  # Using binomial: 1 - P(X = 0) where X ~ Binom(n, pr_detection)
  identity_probs <- sapply(1:3, function(i) {
    sapply(num_cousins[, i], function(n) {1 - pbinom(0, n, pr_detection[i])})
  })
  colnames(identity_probs) <- c("first_cousin", "second_cousin", "third_cousin")
  
  # Return average probability across all simulations
  return(colMeans(identity_probs))
}

#' Simulate Relative Counts Under ZINB Fertility Model
#'
#' Generates simulated counts of relatives at various degrees based on a
#' zero-inflated negative binomial (ZINB) fertility distribution.
#'
#' @param mu Mean parameter of the negative binomial component
#' @param size Dispersion parameter of the negative binomial component
#' @param pi Zero-inflation probability (proportion with no children)
#' @param n_sim Number of individuals to simulate (default: 10000)
#' @param max_kids Maximum number of children to cap at (default: 20)
#' @return Data frame with simulated counts for various relationship types
simulate_relatives_zinb <- function(mu, size, pi,
                                    n_sim = 10000,
                                    max_kids = 20) {
  # ===== Helper Functions =====
  
  # Draw fertility (number of children) from ZINB distribution
  draw_fert <- function(n) {
    # With probability pi, return 0; otherwise draw from negative binomial
    ifelse(runif(n) < pi, 0L,
           as.integer(pmin(rnbinom(n, size, mu = mu), max_kids)))
  }

  # Draw sibling counts from size-biased negative binomial
  draw_sibs <- function(n) {
    as.integer(pmin(rnbinom(n, size + 1, mu = mu * (size + 1) / size), max_kids - 1))
  }

  # ===== Direct Relatives (Fixed or Simple) =====
  
  # Parents and grandparents are deterministic
  parents <- 2L
  grandparents <- 4L

  # Siblings: drawn from size-biased distribution
  siblings <- draw_sibs(n_sim)

  # Children: drawn from fertility distribution
  children <- draw_fert(n_sim)

  # ===== First Degree Relatives =====
  # First degree = parents + siblings + children
  first_degree <- siblings + children + parents

  # ===== Second Degree Relatives =====
  
  # Aunts & uncles: siblings of parents
  # For each parent, draw their number of siblings
  aunts_uncles <- sapply(rep(parents, n_sim), function(n_parents) {
    if (n_parents == 0L) 0L else sum(draw_sibs(n_parents))
  })

  # Nephews & nieces: children of siblings
  # For each sibling, draw their number of children
  nephews_nieces <- sapply(siblings, function(n_sib) {
    if (n_sib == 0L) 0L else sum(draw_fert(n_sib))
  })

  # Second degree = grandparents + aunts/uncles + nephews/nieces
  second_degree <- aunts_uncles + nephews_nieces + grandparents

  # ===== First Cousins =====
  # First cousins: children of aunts/uncles
  # For each aunt/uncle, draw their number of children
  first_cousins <- sapply(aunts_uncles, function(n_au) {
    if (n_au == 0L) 0L else sum(draw_fert(n_au))
  })

  # ===== Second Cousins =====
  # Second cousins: grandchildren of great-aunts/uncles
  # Path: grandparents -> their siblings -> children -> children
  
  # Great-aunts/uncles: siblings of grandparents
  great_aunts_uncles <- sapply(rep(grandparents, n_sim), function(n_grandparents) {
    if (n_grandparents == 0L) 0L else sum(draw_sibs(n_grandparents))
  })
  
  # Children of great-aunts/uncles
  great_aunts_uncles_children <- sapply(great_aunts_uncles, function(n_gau) {
    if (n_gau == 0L) 0L else sum(draw_fert(n_gau))
  })
  
  # Second cousins: children of great-aunts/uncles' children
  second_cousins <- sapply(great_aunts_uncles_children, function(n_gauc) {
    if (n_gauc == 0L) 0L else sum(draw_fert(n_gauc))
  })

  # ===== Third Cousins =====
  # Third cousins: great-grandchildren of great-great-aunts/uncles
  # Path: great-grandparents -> their siblings -> children -> children -> children
  
  # Great-great-aunts/uncles: siblings of great-grandparents (8 total)
  great_great_aunts_uncles <- sapply(rep(grandparents * 2L, n_sim), function(n_ggp) {
    if (n_ggp == 0L) 0L else sum(draw_sibs(n_ggp))
  })
  
  # Children of great-great-aunts/uncles
  great_great_aunts_uncles_children <- sapply(great_great_aunts_uncles, function(n_ggau) {
    if (n_ggau == 0L) 0L else sum(draw_fert(n_ggau))
  })
  
  # Grandchildren of great-great-aunts/uncles
  great_great_aunts_uncles_grandchildren <- sapply(great_great_aunts_uncles_children, function(n_ggauc) {
    if (n_ggauc == 0L) 0L else sum(draw_fert(n_ggauc))
  })
  
  # Third cousins: great-grandchildren of great-great-aunts/uncles
  third_cousins <- sapply(great_great_aunts_uncles_grandchildren, function(n_ggaucg) {
    if (n_ggaucg == 0L) 0L else sum(draw_fert(n_ggaucg))
  })

  # ===== Compile Results =====
  # Return data frame with all simulated relative counts
  relative_counts <- data.frame(num_parents = rep(parents, n_sim),
                                num_siblings = siblings,
                                num_children = children,
                                num_first_degree = first_degree,
                                num_grandparents = rep(grandparents, n_sim),
                                num_aunts_uncles = aunts_uncles,
                                num_nephews_nieces = nephews_nieces,
                                num_second_degree = second_degree,
                                num_first_cousins = first_cousins,
                                num_second_cousins = second_cousins,
                                num_third_cousins = third_cousins
  )
  return(relative_counts)
}


#' Plot Probability Mass Function (PMF)
#'
#' Creates a horizontal bar chart showing the probability mass function
#' for a specified count variable from simulated relative counts.
#'
#' @param relative_counts Data frame containing simulated relative counts
#' @param count_var Name of the column to plot (e.g., "num_siblings")
#' @param category_label Label for the category (e.g., "Siblings")
#' @param cap Optional cap value for grouping larger counts (default: NULL)
#' @return ggplot object showing the PMF as a horizontal bar chart
plot_pmf <- function(relative_counts, count_var = "Siblings", category_label = "Siblings", cap = NULL) {
  # Calculate empirical probability mass function
  pmf_df <- relative_counts %>%
    count(!!rlang::sym(count_var) ) %>%
    mutate(p = n / sum(n))  # Convert counts to probabilities  # Convert counts to probabilities
  
  # Get sorted categories
  cats <- sort(unique(pmf_df[[count_var]]))
  
  # Apply cap if specified (groups all values > cap together)
  if (!is.null(cap)) {
    pmf_df[[count_var]] <- ifelse(pmf_df[[count_var]] > cap, cap, pmf_df[[count_var]])
    cats <- sort(unique(pmf_df[[count_var]]))
  }
  
  # Create color palette using viridis
  n_cat <- length(cats)
  color_palette <- viridis(n_cat, option = "D")
  
  # Convert to factor to preserve ordering
  df <- pmf_df %>%
    mutate(count_cat = factor(!!rlang::sym(count_var), levels = cats))
  
  # Create horizontal bar chart
  ggplot(df, aes(x = count_cat, y = p, fill = count_cat)) +
    geom_col() +  # Bar chart
    geom_hline(yintercept = 0, color = "black", size = 0.5) +  # Baseline
    coord_flip() +  # Horizontal orientation
    scale_y_continuous(labels = function(x) abs(x)) +  # Ensure positive labels
    scale_fill_manual(values = color_palette, name = paste(category_label, "Count")) +
    labs(
      title = paste0("PMF of ", category_label, " Counts (Simulated from ZINB Fertility Model)"),
      x = paste("Number of", category_label),
      y = "Probability"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
      panel.grid.minor.y = element_blank()   # Remove minor grid lines
    )
}
