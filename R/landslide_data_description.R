#' Landslide dataset for the Wenchuan earthquake area
#'
#' This dataset is the outcome of a series of processing steps applied to
#' a global coseismic landslide inventory. The landslide inventory was
#' downloaded from the Global Coseismic Landslide Repository coordinated
#' by the United States Geological Survey (USGS) and described in
#' Tanyaş et al. (2017). The repository provides open-access polygonal
#' and point shapefiles of earthquake-induced landslides.
#'
#' For the Wenchuan earthquake case study, polygonal landslide data were
#' used to compute planimetric landslide areas. The study region was then
#' partitioned into Slope Units, consistent with the spatial units used
#' in Lombardo et al. (2021). Each slope unit was assigned landslide counts,
#' landslide areas, and a set of topographic, hydrological, seismic, and
#' geological predictor variables.
#'
#' @format A data frame with one row per slope unit and the following variables:
#' \describe{
#'   \item{presence}{Binary indicator of landslide occurrence (1 = present, 0 = absent).}
#'   \item{area_grid}{Total area of the spatial grid or slope unit.}
#'   \item{area_slide}{Area of landslide material mapped within the slope unit.}
#'   \item{count}{Number of individual landslide events recorded within the slope unit.}
#'   \item{slope_avg}{Mean slope angle within the slope unit.}
#'   \item{slope_stdd}{Standard deviation of slope within the slope unit.}
#'   \item{relief}{Local terrain relief, defined as the elevation difference within the slope unit.}
#'   \item{TWI_avg}{Mean Topographic Wetness Index value.}
#'   \item{TWI_stddev}{Standard deviation of the Topographic Wetness Index.}
#'   \item{VRM_avg}{Mean Vector Ruggedness Measure, representing surface roughness.}
#'   \item{VRM_stddev}{Standard deviation of the Vector Ruggedness Measure.}
#'   \item{planCurv_a}{Mean plan curvature describing horizontal terrain curvature.}
#'   \item{planCurv_s}{Standard deviation of plan curvature.}
#'   \item{pga_avg}{Mean peak ground acceleration representing seismic shaking intensity.}
#'   \item{pga_stddev}{Standard deviation of peak ground acceleration.}
#'   \item{distStream}{Mean distance to the nearest stream or drainage network.}
#'   \item{distStre_s}{Standard deviation of distance to streams.}
#'   \item{POINT_X}{X-coordinate (longitude or easting) of the slope unit centroid.}
#'   \item{POINT_Y}{Y-coordinate (latitude or northing) of the slope unit centroid.}
#'   \item{litho}{Lithological classification indicating the dominant rock or soil type.}
#'   \item{profCurv_a}{Mean profile curvature describing vertical curvature along the slope direction.}
#'   \item{profCurv_s}{Standard deviation of profile curvature.}
#' }
#'
#' @source
#' Global Coseismic Landslide Repository (USGS):
#' https://usgs.maps.arcgis.com/apps/webappviewer/index.html?id=2b6f1e57135f41028ea42ebc6813d967
#'
#' @references
#' Tanyaş, H., van Westen, C.J., Allstadt, K.E., Jessee, A.N., Görüm, T.,
#' Jibson, R.W., Godt, J.W., Sato, H.P., Schmitt, R.G., Marc, O. and Hovius, N. (2017).
#' Presentation and analysis of a worldwide database of earthquake-induced
#' landslide inventories. Journal of Geophysical Research: Earth Surface,
#' 122(10), 1991--2015.
#'
#' Lombardo, L., Tanyas, H., Huser, R., Guzzetti, F. and Castro-Camilo, D. (2021).
#' Landslide size matters: A new data-driven, spatial prototype.
#' Engineering Geology, 293, 106288.
#'
#' @author
#' Lorenzo Lombardo <l.lombardo@utwente.nl>,
#' Rishikesh Yadav <rishikesh@iitmandi.ac.in>
#'
#' @note
#' No external funding supported the data creation. The authors acknowledge
#' the work of Hakan Tanyas for preparing the original landslide inventory.
#'
#' @docType data
#' @name landslide_data_chapter27
NULL
