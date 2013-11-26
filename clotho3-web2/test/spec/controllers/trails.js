'use strict';

describe('Controller: TrailsCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TrailsCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TrailsCtrl = $controller('TrailsCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
