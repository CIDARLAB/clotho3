'use strict';

Application.Primary.controller('HomeCtrl', ['$scope', '$location', 'Collector', function($scope, $location, Collector) {

    $scope.enterClotho = function() {
        $location.path('/trails/bb99191e810c19729de860fe');
    };

    $scope.enterEugene = function() {
        $location.path('/trails/bb02191e810c19729de860aa');
    };

    $scope.enterRaven = function() {
        $location.path('/trails/bb02191e810c19729de860bb');
    };

    //About page
    $scope.joinMailingList = function(from) {
        Collector.silentAddModel('email_' + Date.now(), from);
        $scope.emailToAdd = null;
    }


}]);