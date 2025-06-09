

function change_year(dt::DateTime, new_year::Int)
    return DateTime(new_year, month(dt), day(dt), hour(dt), minute(dt), second(dt))
end

function change_month(dt::DateTime, new_month::Int)
    return DateTime(year(dt), new_month, day(dt), hour(dt), minute(dt), second(dt))
end