function pln = rmatompln(pln)
%RMATOMPLN Remove atomic planes
    for ipln = numel(pln.plane):-1:1
        if numel(pln.plane(ipln).index) <= 1
            pln.plane(ipln) = [];
        end
    end
end

