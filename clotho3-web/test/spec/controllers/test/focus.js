'use strict';

describe('Controller: TestFocusCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestFocusCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestFocusCtrl = $controller('TestFocusCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
