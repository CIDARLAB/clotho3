'use strict';

describe('Controller: TestSchemaviewCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestSchemaviewCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestSchemaviewCtrl = $controller('TestSchemaviewCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
