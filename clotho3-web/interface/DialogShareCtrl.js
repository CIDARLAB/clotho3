'use strict';

Application.Extensions.controller('DialogShareCtrl', ['$scope', 'dialog', '$location', function($scope, dialog, $location){
    $scope.close = function(result){
        dialog.close(result);
    };

    $scope.social = [
        {
            "name" : "facebook",
            "prefix" : "http://www.facebook.com/sharer.php?u="
        },
        {
            "name" : "google",
            "prefix" : "https://plus.google.com/share?url="
        },
        {
            "name" : "twitter",
            "prefix" : "http://twitter.com/share?url="
        },
        {
            "name" : "linkedin",
            "prefix" : "http://www.linkedin.com/shareArticle?mini=true&url="
        },
        {
            "name" : "digg",
            "prefix" : "http://www.digg.com/submit?url="
        },
        {
            "name" : "reddit",
            "prefix" : "http://reddit.com/submit?url="
        },
        {
            "name" : "email",
            "prefix" : "mailto:?Body="
        }
    ];

    $scope.share = function (site) {
        var url = site.prefix + $location.absUrl();

        $scope.close();

        //join location and site url
        window.open(url, (site.name == 'email' ? '_self' : "_blank") );
    }

}]);