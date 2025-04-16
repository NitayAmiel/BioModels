function label = getModelLabel(version)
%GETMODELLABEL Returns a descriptive label for the given model version.

    switch version
        case 1
            label = 'Baseline model';
        case 2
            label = 'PTM addition';
        case 3
            label = 'Binding affinity reduction';
        otherwise
            label = 'Unknown model';
    end

end
