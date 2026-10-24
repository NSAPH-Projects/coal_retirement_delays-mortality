# Delayed Coal Plant Retirements in the Era of Data Centers Carry a Measurable Mortality Cost

This repository contains the analysis code supporting:

Suhwan Bong, Gianluca Guidi, Michael Cork, Lucas Henneman, and Francesca Dominici. “Delayed Coal Plant Retirements in the Era of Data Centers Carry a Measurable Mortality Cost.” (2026).

## Overview

We compared generator retirement dates reported in the January 2024 and January 2026 EIA-860M inventories. We then estimated plant-specific coal PM2.5 concentrations and mortality under the originally announced and subsequently revised retirement scenarios for 2026–2035.

The study also examines the overlap between the affected plants and balancing authorities containing hyperscale data-center electricity demand.


## Data sources and availability

The repository does not redistribute third-party source data or CMS-restricted Medicare data. The input datasets can be obtained or reconstructed from the sources below.

### EIA generator inventories and retirement dates

Generator-level retirement dates, operating status, plant identifiers, and generator identifiers were obtained from the U.S. Energy Information Administration Preliminary Monthly Electric Generator Inventory (Form
EIA-860M):

https://www.eia.gov/electricity/data/eia860m/

The analysis used the archived January 2024 and January 2026 inventories. The file `coal_plant_retirement_scenarios.csv` was constructed by matching the January 2024 and January 2026 inventories using EIA plant and generator identifiers and applying the retirement-scenario definitions implemented in the analysis code.

### SO2 emissions

Annual plant-level and unit-level SO2 emissions were obtained from the U.S. Environmental Protection Agency Clean Air Markets Program Data:

https://campd.epa.gov/data/custom-data-download

The same records can also be accessed through the CAM API:

https://www.epa.gov/power-sector/cam-api-portal

The file `annual_emission_plant.csv` was constructed from the downloaded CAMPD records.

### Coal PM2.5 source-impact estimates

Historical coal power-plant PM2.5 source-impact estimates were obtained from the long-term coal PM2.5 dataset developed by Henneman et al.:

https://osf.io/8gdau/

Supporting import and processing code is available from:

- [coal_unit_PM25](https://github.com/lhenneman/coal_unit_PM25)
- [coal_exposure_longterm](https://github.com/lhenneman/coal_exposure_longterm)

### ZCTA boundaries and demographic characteristics

2020 Census ZIP Code Tabulation Area boundaries were obtained from the U.S. Census Bureau TIGER/Line files:

https://www.census.gov/geographies/mapping-files/2020/geo/tiger-line-file.html

ZCTA-level demographic characteristics were obtained from the 2020 Census Demographic and Housing Characteristics File:

https://www.census.gov/data/tables/2023/dec/2020-census-dhc.html

### Medicare enrollment and mortality

Medicare enrollment, death counts, and mortality rates were constructed from Centers for Medicare & Medicaid Services research-identifiable administrative data obtained under a data-use agreement.

Information about requesting CMS research-identifiable files is available through ResDAC:

https://resdac.org/research-identifiable-files-rif-requests

The beneficiary-level data and restricted ZCTA-level derivatives cannot be redistributed through this repository. Qualified researchers may request the underlying data through CMS/ResDAC and obtain the required institutional and
CMS approvals.

National Medicare enrollment projections were obtained from the 2025 Annual Report of the Boards of Trustees of the Federal Hospital Insurance and Federal Supplementary Medical Insurance Trust Funds and its supplementary tables:

https://www.cms.gov/data-research/statistics-trends-and-reports/trustees-report-trust-funds

### Balancing-authority assignments

Plant-level balancing-authority assignments were obtained from EPA eGRID2023 Revision 2:

https://www.epa.gov/egrid/detailed-data

EIA plant identifiers were matched to the eGRID `ORISPL` field, and balancing-authority assignments were taken from the corresponding eGRID plant records.

### Hyperscale data-center electricity demand

Balancing-authority-level hyperscale data-center electricity demand was obtained from the analytical framework and accompanying materials of Guidi et al.:

Code and data repository: https://github.com/gianguidi/hyperscale-emissions

The central analysis used the balancing-authority-level HDC demand estimates under the central facility-load scenario.