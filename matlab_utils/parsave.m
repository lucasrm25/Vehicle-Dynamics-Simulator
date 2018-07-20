function parsave(name, var, field)

vars.(field) = var; %#ok<STRNU>
save(name, '-struct', 'vars');

% save(name,'cvs_guma')
