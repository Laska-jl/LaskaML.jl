

"""
    ab_generalized(V, p)

The generalized form af α and β.


"""
function ab_generalized(V, p)
    A, B, C, D, F, H = p
    (A + B * V) / (C + H * exp((V + D) / F))
end
