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
#' Gogé, F., Thuriès, L., Fouad, Y., Damay, N., Davrieux, F., Moussard, G. & Le, R., Caroline & Trupin-Maudemain, Sèverine and Valè, Matthieu & Morvan, T. (2021).
#' Dataset of Chemical and Near-Infrared Spectroscopy Measurements of Fresh and Dried Poultry and Cattle Manure.
#' \emph{Data in Brief}, 34, 106647. \doi{10.1016/j.dib.2020.106647}
#'
#' @examples
#' data(Manure)
"Manure"

#' Sugar Production Data
#'
#' Fluorescence spectroscopy measurements of sugar samples together with
#' laboratory-determined quality parameters (ash content and color).
#'
#' @format A list with two components:
#' \describe{
#'   \item{\code{Y}}{A data frame with 268 observations on the following 2 variables:
#'     \itemize{
#'       \item{\code{Color}}: Absorption at 420 nm of a pH 7 solution,
#'       used as a measure of sugar discoloration (lower values indicate higher quality).
#'       \item{\code{Ash}}: Conductivity-based measure of inorganic
#'       impurities in the refined sugar (percentage).
#'     }}
#'
#'   \item{\code{Fluorescence}}{A 3-way numeric array of dimensions
#'   268 samples x 571 emission wavelengths x 7 excitation wavelengths, where:
#'     \itemize{
#'       \item{Emission wavelengths}: 275–560 nm measured in 0.5 nm intervals.
#'       \item{Excitation wavelengths}: 230, 240, 255, 290, 305, 325, 340 nm.
#'     }
#'     Each entry represents the emission intensity of a sugar solution measured
#'     spectrofluorometrically.}
#' }
#'
#'
#' @references Bro, R. (1999). Exploratory Study of Sugar Production Using Fluorescence Spectroscopy and Multi-way Analysis.
#' \emph{Chemometrics and Intelligent Laboratory Systems}, 46, 133–147.
#' \doi{10.1016/S0169-7439(98)00181-6}
#'
#' @examples
#' data(Sugar)
"Sugar"
