Application.Extensions.controller('clothoIntro_searchbarIntroCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Clotho', function($scope, $focus, $timeout, $dialog, Clotho) {

    $dialog.messageBox('The Clotho Search Bar', 'Clotho\'s Search Bar is a powerful tool, able to do far more than search, and useable anywhere on the site. For example, you could work with data, run commands, or . And yes, you can still search. ', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
    .then(function() {
        $focus.addBackdrop();
        return $focus.highlightElement($('#clothoSearchbar'));
    })
    .then(function() {
        return $focus.removeBackdrop();
    })

}]);