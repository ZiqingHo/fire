#' England Regional Housing Prices
#'
#' A comprehensive dataset of monthly housing prices across England regions.
#'
#' @format A list containing:
#' \describe{
#'   \item{\code{X}}{A 3D array of actual monthly prices in £ (132 months × 4 house types × 9 regions)}
#'   \item{\code{y}}{London's next month's average prices across all property types (£)}
#' }
#'
#' @details
#' \code{X}:
#' \itemize{
#'   \item \emph{Time Period}: January 2000 - December 2010 (132 months)
#'   \item \emph{House Types}:
#'     \itemize{
#'       \item Detached
#'       \item Semi-Detached
#'       \item Terraced
#'       \item Flats
#'     }
#'   \item \emph{Regions}:
#'     \itemize{
#'       \item London
#'       \item East of England
#'       \item South East
#'       \item West Midlands
#'       \item South West
#'       \item East Midlands
#'       \item North West
#'       \item Yorkshire & Humber
#'       \item North East
#'     }
#' }
#'
#' \code{y}:
#' \itemize{
#'   \item Month: From February 2000 to January 2011
#'   \item Price: Next month's average prices across all property types in London
#' }
#'
#' @source HM Land Registry Open Data (2024). UK House Price Index.
#' \url{https://landregistry.data.gov.uk/#ukhpi}.
#' Downloaded on 2024-12-24 and processed by package author.
#'
#' @examples
#' data(Housing)
#'
#' # Get all detached home prices in London
#' london_detached <- Housing$X[, "Detached", "London"]
"Housing"

#' French Manure Data
#'
#' Combined near-infrared spectra and chemical composition measurements of cattle and poultry manure samples.
#'
#' @format A list containing two components:
#' \describe{
#'   \item{\code{absorp}}{A matrix of NIR absorbance spectra (332 samples × 700 wavelengths). Wavelength range: 1100-2500 nm (2 nm resolution)}
#'   \item{\code{y}}{A matrix of chemical measurements (332 samples × 3 variables)}
#'   \itemize{
#'   \item \code{DM}: Dry matter content (\% of wet weight)
#'   \item \code{NH4}: Total ammonium nitrogen (\% of wet weight)
#'   \item \code{N}: Total nitrogen (\% of wet weight)
#' }
#' }
#'
#' @references
#' Gogé, F., Thuriès, L., Fouad, Y., Damay, N., Davrieux, F., Moussard, G., ... & Morvan, T. (2021).
#' Dataset of chemical and near-infrared spectroscopy measurements of fresh and dried poultry and cattle manure.
#' \emph{Data in Brief}, 34, 106647. \doi{10.1016/j.dib.2020.106647}
#'
#' @examples
#' data(Manure)
"Manure"

#' PeMS Highway Traffic Data
#'
#' Hourly traffic speed measurements from California highway sensors.
#'
#' @format A 3-way tensor (array) with dimensions:
#' \describe{
#'   \item{Hours}{24 hours of daily measurements}
#'   \item{Days}{44 days of observations (weekdays of May-June 2012)}
#'   \item{Sensors}{288 sensor stations in District 7}
#' }
#'
#' @details
#' Processed from the original PeMSD7 dataset containing raw 5-minutes measurements aggregated to hourly intervals
#'
#' @source
#' Downloaded on 2025-07-19 and processed by package author.
#'
#' @references
#' B. Yu, H. Yin, and Z. Zhu, "Spatio-temporal graph convolutional networks:
#' A deep learning framework for traffic forecasting," in Proc. IJCAI, 2018, pp. 3634–3640.
#'
#' @examples
#' data(PeMS)
"PeMS"


