import NProgress from 'nprogress';

$( document ).ajaxStart(function() {
      NProgress.start();
});


$( document ).ajaxStop(function() {
     NProgress.done();
});
