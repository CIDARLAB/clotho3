'use strict';

Application.Extensions.controller('TestDialogController', function($scope, dialog){
    $scope.close = function(result){
        dialog.close(result);
    };

    $scope.filterTest = "Here is some Text";
});