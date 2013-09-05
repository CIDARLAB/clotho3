'use strict';

Application.Extensions.controller('clothoIntro_scriptingIntroCtrl', ['$scope', '$focus', '$timeout', function($scope, $focus, $timeout) {

    $scope.showMeHow = function() {
        var searchInput = $('#searchBarInput'),
            searchSubmit = $('#searchBarSubmit'),
            oldinputZ = searchInput.css("z-index"),
            oldsubmitZ = searchSubmit.css("z-index"),
            maxZ = $focus.maxZ;

        //setup
        searchInput.val('');
        searchSubmit.css('position', 'relative');

        //add backdrop, highlight input
        searchInput.css("z-index", maxZ + 2);
        $focus.addBackdrop(maxZ+1);

        $timeout(function() {
            //start typing
            searchInput.val("var result = clotho.run('revcomp', ['acagtgcca'])");
        }, 500)
        .then(function() {
            return $timeout(function() {
                //fade out search
                searchInput.css("z-index", oldinputZ);
            }, 500)
        })
        .then(function() {
            return $timeout(function() {
                //focus submit
                searchSubmit.css("z-index", maxZ + 2);
            }, 100)
        })
        .then(function() {
            return $timeout(function() {
                //submit

                //fixme
                searchSubmit.click();

                $focus.removeBackdrop();
                searchSubmit.css("z-index", oldsubmitZ)
            }, 500)
        });
    }
}]);