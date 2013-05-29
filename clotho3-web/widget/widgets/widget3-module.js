'use strict';

angular.module('widgetApp3', ['clotho.core'])
    //todo - rewrite so still function but not reliant on route
    .config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            otherwise({
                templateUrl:'widget/widgets/widget3-partial.html'
            })
    }])
    .run(['$rootScope', 'Clotho', function($rootScope, Clotho) {
        $rootScope.Clotho = Clotho;
    }])
    .controller('widget3Ctrl', ['$scope', 'Clotho', function($scope, Clotho) {
        $scope.greeting = "this is a controller string";

        $scope.menuItems = Clotho.get('menu_items');
    }]);

