'use strict';

describe('Controller: TestContstructiontrailCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestContstructiontrailCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestContstructiontrailCtrl = $controller('TestContstructiontrailCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
