'use strict';

describe('Controller: TestConstructionCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestConstructionCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestConstructionCtrl = $controller('TestConstructionCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
