"""
contains the information of regions in the U.S.
"""


const g_statesUS48 = ["New Jersey","Rhode Island","Massachusetts","Connecticut","Maryland","New York","Delaware","Florida","Ohio","Pennsylvania","Illinois","California","Virginia","Michigan","Indiana","North Carolina","Georgia","Tennessee","New Hampshire","South Carolina","Louisiana","Kentucky","Wisconsin","Washington","Alabama","Missouri","Texas","West Virginia","Vermont","Minnesota","Mississippi","Iowa","Arkansas","Oklahoma","Arizona","Colorado","Maine","Oregon","Kansas","Utah","Nebraska","Nevada","Idaho","New Mexico","South Dakota","North Dakota","Montana","Wyoming",
] # "Alaska",  "Hawaii",



g_us_state_abbr_dict = Dict(
    "Alabama" => "AL",
    "Alaska" => "AK",
    "Arizona" => "AZ",
    "Arkansas" => "AR",
    "California" => "CA",
    "Colorado" => "CO",
    "Connecticut" => "CT",
    "Delaware" => "DE",
    "District of Columbia" => "DC",
    "Florida" => "FL",
    "Georgia" => "GA",
    "Hawaii" => "HI",
    "Idaho" => "ID",
    "Illinois" => "IL",
    "Indiana" => "IN",
    "Iowa" => "IA",
    "Kansas" => "KS",
    "Kentucky" => "KY",
    "Louisiana" => "LA",
    "Maine" => "ME",
    "Maryland" => "MD",
    "Massachusetts" => "MA",
    "Michigan" => "MI",
    "Minnesota" => "MN",
    "Mississippi" => "MS",
    "Missouri" => "MO",
    "Montana" => "MT",
    "Nebraska" => "NE",
    "Nevada" => "NV",
    "New Hampshire" => "NH",
    "New Jersey" => "NJ",
    "New Mexico" => "NM",
    "New York" => "NY",
    "North Carolina" => "NC",
    "North Dakota" => "ND",
    "Ohio" => "OH",
    "Oklahoma" => "OK",
    "Oregon" => "OR",
    "Pennsylvania" => "PA",
    "Puerto Rico" => "PR",
    "Rhode Island" => "RI",
    "South Carolina" => "SC",
    "South Dakota" => "SD",
    "Tennessee" => "TN",
    "Texas" => "TX",
    "Utah" => "UT",
    "Vermont" => "VT",
    "Virginia" => "VA",
    "Washington" => "WA",
    "West Virginia" => "WV",
    "Wisconsin" => "WI",
    "Wyoming" => "WY"
)

g_us_abbr48_vec = [g_us_state_abbr_dict[state] for state in g_statesUS48]


g_us_abbr_to_state_dict = Dict(values(g_us_state_abbr_dict) .=> keys(g_us_state_abbr_dict))

