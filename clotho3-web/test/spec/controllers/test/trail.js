'use strict';

describe('Controller: TestTrailCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestTrailCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestTrailCtrl = $controller('TestTrailCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
