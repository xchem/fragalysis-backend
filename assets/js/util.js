// return querystring parameter by name
function getQueryStringParameterByName(name) {
    const href = window.location.href;
    name = name.replace(/[\[\]]/g, "\\$&");
    const regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)");
    const results = regex.exec(href);
    //if (!results) return null;
    if (!results) return '';  // return empty string so we don't have to check for null
    if (!results[2]) return '';
    return decodeURIComponent(results[2].replace(/\+/g, " "));
}

// inspect any obj for debugging - TODO add recursion
function inspect(obj) {
    console.log('util: inspect --------------------------');
    if (Array.isArray(obj)) { // array is obj too so expicitly check for it before checking for object
        console.log('obj is array');
        obj.forEach((item) => console.log(item));
    } else if (obj instanceof Object){
        for (var prop in obj) {
            console.log(prop + '=' + obj[prop]);
        }
    } else {
        console.log(typeof obj);
    }
}

export { getQueryStringParameterByName, inspect }
