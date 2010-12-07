function post = get_posted_sel(obj)

% GET_POSTED_SEL get_posted_sel(pdf4D) returns posted selections

%test if msi was installed
if ~msi_installed
    error('MSI was not installed')
end

%get posted selection from msi
[stat, post] = unix('get_posted_sel');

%in case get_posted_sel returns error
if stat~=0
    error('Can not get posted selection')
end
