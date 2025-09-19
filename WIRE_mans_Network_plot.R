#install.packages("rlang")
#packageVersion("rlang")
#install.packages("magrittr", type = "binary")
#install.packages("ergm", dependencies = TRUE, type = "binary")
#install.packages("statnet",dependencies = TRUE, type = "binary")
#sessionInfo()


library(ergm)
library(statnet)
library(sna)
library(tergm)
library(ergm.count)
library(statnet.common)
library(tsna)
library(network)
library(networkDynamic)


################## Dynamic Network #################
####################################################

library(igraph)

# Time 1: baseline karate network
g1 <- make_graph("Zachary")

# Assign names if missing
if (is.null(V(g1)$name)) {
  V(g1)$name <- as.character(V(g1))  # names = "1","2",...
}

# Time 2: modify snapshot
set.seed(42)
rem_ids <- sample(E(g1), 6)
g2 <- delete_edges(g1, rem_ids)

to_add <- matrix(c(
  1, 5,
  9, 34,
  14, 20,
  24, 33,
  10, 28,
  3, 17
), ncol = 2, byrow = TRUE)

need_add <- !apply(to_add, 1, function(p) are_adjacent(g2, p[1], p[2]))
to_add2  <- to_add[need_add, , drop = FALSE]
if (nrow(to_add2) > 0) g2 <- add_edges(g2, t(to_add2))

# Ensure g2 also has names
if (is.null(V(g2)$name)) {
  V(g2)$name <- as.character(V(g2))
}


# --- Shared layout for visual comparability ---
set.seed(99)
lay <- layout_with_fr(g1)

# --- Side-by-side plot ---
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

plot(g1,
     layout = lay,
     vertex.size = 18,
     vertex.label = V(g1)$name,         # show vertex names
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.color = "skyblue",
     edge.color = "grey40",
     main = "Time 1")

plot(g2,
     layout = lay,
     vertex.size = 18,
     vertex.label = V(g2)$name,         # show vertex names
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.color = "lightgreen",
     edge.color = "grey40",
     main = "Time 2")

par(mfrow = c(1, 1))

# Now union works
gU <- igraph::union(g1, g2, byname = TRUE)

# Helper to form an undirected edge key "u--v" with sorted endpoints
ekey <- function(g) {
  apply(ends(g, E(g)), 1, function(x) paste(sort(x), collapse = "--"))
}

e1 <- ekey(g1)
e2 <- ekey(g2)
eU <- ekey(gU)

# Classify edges by change status
persist   <- intersect(e1, e2)
t1_only   <- setdiff(e1, e2)
t2_only   <- setdiff(e2, e1)

edge_col <- ifelse(eU %in% persist, "steelblue3",
                   ifelse(eU %in% t1_only, "skyblue", "lightgreen"))
edge_lty <- ifelse(eU %in% persist, 1,
                   ifelse(eU %in% t1_only, 2, 1))
edge_lwd <- ifelse(eU %in% persist, 1.6,
                   ifelse(eU %in% t1_only, 2.2, 2.2))

# Single layout for fair comparison
set.seed(99)
lay <- layout_with_fr(gU)

# Optional: size by degree in the union (stable visual emphasis)
vdeg <- degree(gU)

# Plot
plot(
  gU,
  layout       = lay,
  vertex.size  = scales::rescale(vdeg, to = c(10, 18)),
  vertex.label = V(gU)$name,
  vertex.label.cex = 0.8,
  vertex.frame.color = "black",
  vertex.color = "steelblue3",
  edge.color   = edge_col,
  edge.lty     = edge_lty,
  edge.width   = edge_lwd,
  main = "Dynamic Network Representation"
)

legend(
  "topright",
  legend = c("Persists", "Only at Time 1", "Only at Time 2"),
  lty    = c(1, 2, 1),
  col    = c("steelblue3", "skyblue", "lightgreen"),
  lwd    = c(1.6, 2.2, 2.2),
  bty    = "n"
)


################## Global and Local Representation ###########
##############################################################

## --- Data ---
g <- make_graph("Zachary")
V(g)$name <- as.character(V(g))  # ensure names 1..34 as character

## --- Global styling: communities & degree ---
set.seed(1)
lay <- layout_with_fr(g)                 # single shared layout
deg <- degree(g)
com <- cluster_fast_greedy(as_undirected(g))
pal <- grDevices::hcl.colors(length(unique(membership(com))), "Set3")
vcol <- pal[membership(com)]

## --- Choose a focal node and ego radius (local view) ---
focal <- "34"         # change to "34" or any vertex name
radius <- 2          # 1-hop ego; set to 2 for a larger “local” region

ego_vids <- ego(g, order = radius, nodes = focal, mode = "all")[[1]]
g_local  <- induced_subgraph(g, ego_vids)

## Keep the SAME coordinates for the local nodes (subset the global layout)
# Map row order: vertex sequence of g_local vs g
idx <- match(V(g_local)$name, V(g)$name)
lay_local <- lay[idx, , drop = FALSE]

## --- Plot: side by side ---
op <- par(mfrow = c(1, 2), mar = c(1,1,3,1))

# Global
plot(g,
     layout = lay,
     vertex.size = scales::rescale(deg, to = c(8, 18)),
     vertex.label = V(g)$name,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.color = vcol,
     vertex.frame.color = NA,
     edge.color = "grey60",
     main = "Global network (communities, size = degree)")

# Local (ego)
hl <- ifelse(V(g)$name %in% V(g_local)$name, "gold", "grey85")  # highlight ego nodes
plot(g_local,
     layout = lay_local,
     vertex.size = scales::rescale(deg[idx], to = c(10, 20)),
     vertex.label = V(g_local)$name,
     vertex.label.cex = 0.9,
     vertex.label.color = "black",
     vertex.color = "skyblue",
     vertex.frame.color = "grey30",
     edge.color = "steelblue3",
     main = sprintf("Local structure: ego(node=%s, radius=%d)", focal, radius))

par(op)







