'use strict';

Application.Extensions.controller('clothoIntro_scriptingIntroCtrl', ['$scope', '$focus', '$timeout', function($scope, $focus, $timeout) {

    $scope.showMeHow = function() {
        var searchInput = $('#searchBarInput'),
            searchSubmit = $('#searchBarSubmit'),
            oldinputZ = searchInput.css("z-index"),
            oldsubmitZ = searchSubmit.css("z-index"),
            maxZ = $focus.maxZ;


        //add backdrop, highlight input
        $focus.setZ(maxZ + 2, searchInput)
        .then(function() {
            return $focus.addBackdrop(maxZ+1)
        })
        .then(function() {
            return $focus.typeOut(searchInput,
                "clotho.run('revcomp', ['acagtgcca'])", 'display.query')
        })
        .then(function() {
            return $timeout(function() {
                //fade out search
                $focus.setZ(oldinputZ, searchInput);
            }, 500)
        })
        .then(function() {
            return $timeout(function() {
                //focus submit
                $focus.setZ(maxZ + 2, searchSubmit);
            }, 100)
        })
        .then(function() {
            return $timeout(function() {
                //submit
                searchSubmit.click();

                $focus.removeBackdrop();
                $focus.setZ(oldsubmitZ, searchSubmit);
            }, 500)
        });
    }
}]);