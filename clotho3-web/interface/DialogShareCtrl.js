'use strict';

Application.Extensions.controller('DialogShareCtrl', function($scope, $dialog){
    $scope.close = function(result){
        $dialog.close(result);
    };

    //todo - implement this and share()
    $scope.social = [
        {
            "name" : "facebook",
            "url" : ""
        },
        {
            "name" : "google",
            "url" : ""
        },
        {
            "name" : "twitter",
            "url" : ""
        }
    ];

    $scope.share = function (site) {
        //join location and site url
        alert('you would share this page on ' + site.name + '....');
    }

});